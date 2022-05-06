#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <filesystem>
#include <vector>
#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <thread>
#include <boost/format.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <mongocxx/instance.hpp>
#include <mongocxx/pool.hpp>
#include <mongocxx/client.hpp>
#include <bsoncxx/json.hpp>
#include "utility.hpp"
#include "safe_counter.hpp"
#include "io_service_pool.hpp"
#include "compound_database.hpp"
using namespace std;
using namespace std::chrono;
using namespace std::filesystem;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::MorganFingerprints;
using namespace RDGeom;
using namespace RDNumeric::Alignments;
using namespace MolTransforms;
using namespace RDKit::Descriptors;
using namespace mongocxx;
using bsoncxx::builder::basic::kvp;

template<typename T>
auto dist2(const T& p0, const T& p1)
{
	const auto d0 = p0[0] - p1[0];
	const auto d1 = p0[1] - p1[1];
	const auto d2 = p0[2] - p1[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

array<Point3D, 4> calcRefPoints(const ROMol& mol, const vector<int>& heavyAtoms)
{
	const auto num_points = heavyAtoms.size();
	assert(num_points == mol.getNumHeavyAtoms());
	const auto& conf = mol.getConformer();
	array<Point3D, 4> refPoints;
	for (auto& ref : refPoints)
	{
		assert(ref[0] == 0);
		assert(ref[1] == 0);
		assert(ref[2] == 0);
	}
	auto& ctd = refPoints[0];
	auto& cst = refPoints[1];
	auto& fct = refPoints[2];
	auto& ftf = refPoints[3];
	for (const auto i : heavyAtoms)
	{
		const auto& a = conf.getAtomPos(i);
		ctd += a;
	}
	ctd /= static_cast<double>(num_points);
	auto cst_dist = numeric_limits<double>::max();
	auto fct_dist = numeric_limits<double>::lowest();
	auto ftf_dist = numeric_limits<double>::lowest();
	for (const auto i : heavyAtoms)
	{
		const auto& a = conf.getAtomPos(i);
		const auto this_dist = dist2(a, ctd);
		if (this_dist < cst_dist)
		{
			cst = a;
			cst_dist = this_dist;
		}
		if (this_dist > fct_dist)
		{
			fct = a;
			fct_dist = this_dist;
		}
	}
	for (const auto i : heavyAtoms)
	{
		const auto& a = conf.getAtomPos(i);
		const auto this_dist = dist2(a, fct);
		if (this_dist > ftf_dist)
		{
			ftf = a;
			ftf_dist = this_dist;
		}
	}
	return refPoints;
}

int main(int argc, char* argv[])
{
	// Check the required number of command line arguments.
	if (argc != 6)
	{
		cout << "lbvs databases host port user pass" << endl;
		return 0;
	}

	// Fetch command line arguments.
	const path cpdbs_path = argv[1];
	const string host = argv[2];
	const string port = argv[3];
	const string user = argv[4];
	const string pass = argv[5];

	// Initialize constants.
	cout << local_time_string() << "Initializing" << endl;
	const size_t num_usrs = 2;
	const array<string, 2> usr_names{{ "USR", "USRCAT" }};
	constexpr array<size_t, num_usrs> qn{{ 12, 60 }};
	constexpr array<double, num_usrs> qv{{ 1.0 / qn[0], 1.0 / qn[1] }};
	const size_t num_refPoints = 4;
	const size_t num_subsets = 5;
	const array<string, num_subsets> SubsetSMARTS
	{{
		"[!#1]", // heavy
		"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
		"[a]", // aromatic
		"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
		"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
	}};
	const size_t num_hits = 100;

	// Wrap SMARTS strings to RWMol objects.
	array<unique_ptr<ROMol>, num_subsets> SubsetMols;
	for (size_t k = 0; k < num_subsets; ++k)
	{
		SubsetMols[k].reset(reinterpret_cast<ROMol*>(SmartsToMol(SubsetSMARTS[k])));
	}

	// Read compound database directory.
	vector<compound_database> databases;
	for (directory_iterator cpdbs_dir_iter(cpdbs_path), end_dir_iter; cpdbs_dir_iter != end_dir_iter; ++cpdbs_dir_iter)
	{
		// Filter out non directory.
		if (!cpdbs_dir_iter->is_directory()) continue;

		// Create a compound database instance.
		databases.emplace_back(cpdbs_dir_iter->path());
	}
	string db_op_in_array;
	for (size_t i = 0; i < databases.size(); ++i)
	{
		db_op_in_array += (i ? ",\"" : "\"") + databases[i].name + "\"";
	}

	// Connect to mongodb and authenticate user.
	cout << local_time_string() << "Connecting to " << host << ':' << port << " and authenticating " << user << endl;
	const instance inst; // The constructor and destructor initialize and shut down the driver. http://mongocxx.org/api/current/classmongocxx_1_1instance.html
	const uri uri("mongodb://" + host + ":" + port + "/?minPoolSize=0&maxPoolSize=2"); // When connecting to a replica set, it is much more efficient to use a pool as opposed to manually constructing client objects. TODO: Create a credential that will authenticate properly regardless of server version, use a connection string with the user and password directly in the URI and with a parameter specifying the database to authenticate from: "mongodb://user:pass@localhost:27017/?authSource=jstar"
	pool pool(uri);
	const auto client = pool.acquire(); // Return value of acquire() is an instance of entry. An entry is a handle on a client object acquired via the pool.
	auto coll = client->database("jstar").collection("lbvs");
	const auto jobid_filter = bsoncxx::from_json(R"({ "startDate" : { "$exists" : false }, "database": { "$in": [)" + db_op_in_array + R"(] }})");
	const auto jobid_foau_options = options::find_one_and_update().sort(bsoncxx::from_json(R"({ "submitDate" : 1 })")).projection(bsoncxx::from_json(R"({ "_id" : 1, "qryMolSdf": 1, "database": 1, "score": 1 })")); // By default, the original document is returned

	// Initialize variables.
	array<vector<int>, num_subsets> subsets;
	array<vector<double>, num_refPoints> dista;
	alignas(32) array<double, 60> q; // The usrcat feature vector of the query molecule.

	// Initialize an io service pool and create worker threads for later use.
	const size_t num_threads = thread::hardware_concurrency();
	cout << local_time_string() << "Creating an io service pool of " << num_threads << " worker threads" << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;

	// Enter event loop.
	cout << local_time_string() << "Entering event loop" << endl;
	cout.setf(ios::fixed, ios::floatfield);
	bool sleeping = false;
	while (true)
	{
		// Fetch an incompleted job in a first-come-first-served manner.
		if (!sleeping) cout << local_time_string() << "Fetching an incompleted job" << endl;
		const auto startDate = system_clock::now();
		bsoncxx::builder::basic::document jobid_update_builder;
		jobid_update_builder.append(
			kvp("$set", [=](bsoncxx::builder::basic::sub_document set_subdoc) {
				set_subdoc.append(kvp("startDate", bsoncxx::types::b_date(startDate)));
			})
		);

		// http://mongocxx.org/api/current/classmongocxx_1_1collection.html#a6b04b5265e044c08d1ecd1e4185fb699
		const auto jobid_document = coll.find_one_and_update(jobid_filter.view(), jobid_update_builder.extract(), jobid_foau_options);
		if (!jobid_document)
		{
			// No incompleted jobs. Sleep for a while.
			if (!sleeping) cout << local_time_string() << "Sleeping" << endl;
			sleeping = true;
			this_thread::sleep_for(std::chrono::seconds(2));
			continue;
		}
		sleeping = false;
		const auto jobid_view = jobid_document->view();

		// Obtain job properties.
		const auto _id = jobid_view["_id"].get_oid().value;
		cout << local_time_string() << "Executing job " << _id.to_string() << endl;
		const auto qry_mol_sdf = jobid_view["qryMolSdf"].get_utf8().value; // get_utf8().value returns an instance of std::string_view.
		const auto cpdb_name = jobid_view["database"].get_utf8().value;
		const auto score = jobid_view["score"].get_utf8().value;
		assert(score.compare("USR") == 0 || score.compare("USRCAT") == 0);
		const size_t usr0 = score.compare("USRCAT") == 0; // Specify the primary sorting score. 0: USR; 1: USRCAT.
		const auto usr1 = usr0 ^ 1;
		const auto qnu0 = qn[usr0];
		const auto qnu1 = qn[usr1];

		// Obtain a constant reference to the selected database.
		cout << local_time_string() << "Finding the selected compound database" << endl;
		size_t cpdb_index = 0;
		while (cpdb_name.compare(databases[cpdb_index].name) != 0) ++cpdb_index;
		const auto& cpdb = databases[cpdb_index];

		// Read the user-supplied SDF file.
		cout << local_time_string() << "Reading the query file" << endl;
		istringstream qry_mol_sdf_iss(qry_mol_sdf.data()); // data() may return a pointer to a buffer that is not null-terminated. Therefore it is typically a mistake to pass data() to a routine that takes just a const CharT* and expects a null-terminated string.
		SDMolSupplier qry_mol_sup(&qry_mol_sdf_iss, false, true, false, true); // takeOwnership, sanitize, removeHs, strictParsing. Note: setting removeHs=true (which is the default setting) will lead to fewer hydrogen bond acceptors being matched.

		// Initialize vectors to store compounds' primary score and their corresponding conformer.
		vector<double> scores(cpdb.num_compounds); // Primary score of compounds.
		vector<size_t> cnfids(cpdb.num_compounds); // ID of conformer with the best primary score.
		const auto compare = [&](const size_t val0, const size_t val1) // Sort by the primary score.
		{
			return scores[val0] < scores[val1];
		};

		// Initialize the number of chunks and the number of compounds per chunk.
		const auto chunk_size = max(1 + (cpdb.num_compounds - 1) / (num_threads << 2), num_hits);
		const auto num_chunks = 1 + (cpdb.num_compounds - 1) / chunk_size;
		assert(chunk_size * num_chunks >= cpdb.num_compounds);
		assert(chunk_size >= num_hits);
		cout << local_time_string() << "Using " << num_chunks << " chunks and a chunk size of " << chunk_size << endl;
		vector<size_t> scase(cpdb.num_compounds);
		vector<size_t> zcase(num_hits * (num_chunks - 1) + min(num_hits, cpdb.num_compounds - chunk_size * (num_chunks - 1))); // The last chunk might have fewer than num_hits records.
		ifstream conformers_sdf_ifs(cpdb.dpth / "conformers.sdf");

		// Process each of the query compounds sequentially.
		ostringstream hit_mol_sdf_oss;
		const auto num_qry_mols = qry_mol_sup.length(); // num_qry_mols is the number of query molecules submitted by the user. These query molecules are not necessarily all processed, given the limitation of maximum 16MB MongoDB document size of the result.
		unsigned int query_number = 0;
		while (query_number < num_qry_mols)
		{
			cout << local_time_string() << "Parsing query compound " << query_number << endl;
			const unique_ptr<ROMol> qry_mol_ptr(qry_mol_sup.next()); // Calling next() may print "ERROR: Could not sanitize compound on line XXXX" to stderr.
			auto& qryMol = *qry_mol_ptr;

			// Get the number of atoms, including and excluding hydrogens.
			const auto num_atoms = qryMol.getNumAtoms();
			const auto num_heavy_atoms = qryMol.getNumHeavyAtoms();
			assert(num_heavy_atoms);
			cout << local_time_string() << "Found " << num_atoms << " atoms and " << num_heavy_atoms << " heavy atoms" << endl;

			// Calculate Morgan fingerprint.
			cout << local_time_string() << "Calculating Morgan fingerprint" << endl;
			const unique_ptr<SparseIntVect<uint32_t>> qryFp(getFingerprint(qryMol, 2));

			// Classify atoms to pharmacophoric subsets.
			cout << local_time_string() << "Classifying atoms into subsets" << endl;
			for (size_t k = 0; k < num_subsets; ++k)
			{
				vector<vector<pair<int, int>>> matchVect;
				SubstructMatch(qryMol, *SubsetMols[k], matchVect);
				const auto num_matches = matchVect.size();
				auto& subset = subsets[k];
				subset.resize(num_matches);
				for (size_t i = 0; i < num_matches; ++i)
				{
					subset[i] = matchVect[i].front().second;
				}
				cout << local_time_string() << "Found " << num_matches << " atoms for subset " << k << endl;
			}
			const auto& subset0 = subsets.front();
			assert(subset0.size() == num_heavy_atoms);

			// Calculate the four reference points.
			cout << local_time_string() << "Calculating " << num_refPoints << " reference points" << endl;
			const auto qryRefPoints = calcRefPoints(qryMol, subset0);
			const Point3DConstPtrVect qryRefPointv
			{{
				&qryRefPoints[0],
				&qryRefPoints[1],
				&qryRefPoints[2],
				&qryRefPoints[3],
			}};

			// Precalculate the distances of heavy atoms to the reference points, given that subsets[1 to 4] are subsets of subsets[0].
			cout << local_time_string() << "Calculating " << num_heavy_atoms * num_refPoints << " pairwise distances" << endl;
			const auto& qryCnf = qryMol.getConformer();
			for (size_t k = 0; k < num_refPoints; ++k)
			{
				const auto& refPoint = qryRefPoints[k];
				auto& distp = dista[k];
				distp.reserve(num_atoms);
				for (size_t i = 0; i < num_heavy_atoms; ++i)
				{
					distp[subset0[i]] = sqrt(dist2(qryCnf.getAtomPos(subset0[i]), refPoint));
				}
			}

			// Loop over pharmacophoric subsets and reference points.
			cout << local_time_string() << "Calculating " << 3 * num_refPoints * num_subsets << " moments of USRCAT feature" << endl;
			size_t qo = 0;
			for (const auto& subset : subsets)
			{
				const auto n = subset.size();
				for (size_t k = 0; k < num_refPoints; ++k)
				{
					// Load distances from precalculated ones.
					const auto& distp = dista[k];
					vector<double> dists(n);
					for (size_t i = 0; i < n; ++i)
					{
						dists[i] = distp[subset[i]];
					}

					// Compute moments.
					array<double, 3> m{};
					if (n > 2)
					{
						const auto v = 1.0 / n;
						for (size_t i = 0; i < n; ++i)
						{
							const auto d = dists[i];
							m[0] += d;
						}
						m[0] *= v;
						for (size_t i = 0; i < n; ++i)
						{
							const auto d = dists[i] - m[0];
							m[1] += d * d;
						}
						m[1] = sqrt(m[1] * v);
						for (size_t i = 0; i < n; ++i)
						{
							const auto d = dists[i] - m[0];
							m[2] += d * d * d;
						}
						m[2] = cbrt(m[2] * v);
					}
					else if (n == 2)
					{
						m[0] = 0.5 *     (dists[0] + dists[1]);
						m[1] = 0.5 * fabs(dists[0] - dists[1]);
					}
					else if (n == 1)
					{
						m[0] = dists[0];
					}
					for (const auto e : m)
					{
						q[qo++] = e;
					}
				}
			}
			assert(qo == qn.back());

			// Compute USR and USRCAT scores.
			cout << local_time_string() << "Screening " << cpdb.name << " and calculating " << cpdb.num_compounds << " " << usr_names[usr0] << " scores from " << cpdb.num_conformers << " conformers" << endl;
			scores.assign(scores.size(), numeric_limits<double>::max());
			iota(scase.begin(), scase.end(), 0);
			cnt.init(num_chunks);
			for (size_t l = 0; l < num_chunks; ++l)
			{
				io.post([&,l]()
				{
					// Loop over compounds of the current chunk.
					const auto chunk_beg = chunk_size * l;
					const auto chunk_end = min(chunk_beg + chunk_size, cpdb.num_compounds);
					for (size_t k = chunk_beg; k < chunk_end; ++k)
					{
						// Loop over conformers of the current compound and calculate their primary score.
						auto& scorek = scores[k];
						for (size_t j = k << 2; j < (k + 1) << 2; ++j)
						{
							const auto& d = cpdb.usrcat[j];
							double s = 0;
							for (size_t i = 0; i < qnu0; ++i)
							{
								s += abs(q[i] - d[i]);
								if (s >= scorek) break;
							}
							if (s < scorek)
							{
								scorek = s;
								cnfids[k] = j;
							}
						}
					}

					// Sort the scores of compounds of the current chunk.
					sort(scase.begin() + chunk_beg, scase.begin() + chunk_end, compare);

					// Copy the indexes of top hits of the current chunk to a global vector for final sorting.
					copy_n(scase.begin() + chunk_beg, min(num_hits, chunk_end - chunk_beg), zcase.begin() + num_hits * l);

					cnt.increment();
				});
			}
			cnt.wait();

			// Sort the top hits from chunks.
			cout << local_time_string() << "Sorting " << zcase.size() << " hits by " << usr_names[usr0] << " score" << endl;
			sort(zcase.begin(), zcase.end(), compare);

			// Write hit molecules to a string stream for output.
			// This loop can be parallelized, but beware of using independent conformers_sdf_ifs.
			// 1) array<ostringstream, num_hits> hit_mol_sdf_oss_arr; cnt.init(num_hits);
			// 2) for (size_t l = 0; l < num_hits; ++l) { io.post([&,l](){ SDWriter hit_mol_sdf_writer(&hit_mol_sdf_oss_arr[l], false); }); cnt.increment(); }
			// 3) cnt.wait(); ostringstream hit_mol_sdf_oss; for (size_t l = 0; l < num_hits; ++l) { hit_mol_sdf_oss << hit_mol_sdf_oss_arr[l]; } const auto hit_mol_sdf = hit_mol_sdf_oss.str();
			cout << local_time_string() << "Writing hit molecules to a string stream" << endl;
			ostringstream hit_mol_sdf_per_qry_oss;
			SDWriter hit_mol_sdf_writer(&hit_mol_sdf_per_qry_oss, false); // std::ostream*, bool takeOwnership
			for (size_t l = 0; l < num_hits; ++l)
			{
				// Obtain indexes to the hit compound and the hit conformer.
				const auto k = zcase[l];
				const auto j = cnfids[k];
				cout << l << ' ' << k << ' ' << j << endl;

				// Calculate the secondary score of the saved conformer, which has the best primary score.
				const auto& d = cpdb.usrcat[j];
				double s = 0;
				for (size_t i = 0; i < qnu1; ++i)
				{
					s += abs(q[i] - d[i]);
				}
				const auto u0score = 1 / (1 + scores[k] * qv[usr0]); // Primary score of the current compound.
				const auto u1score = 1 / (1 + s * qv[usr1]); // Secondary score of the current compound.

				// Read SDF content of the hit conformer.
				istringstream hit_mol_sdf_iss(cpdb.read_conformer(j, conformers_sdf_ifs));

				// Construct a RDKit ROMol object.
				SDMolSupplier hit_mol_sup(&hit_mol_sdf_iss, false, true, false, true);
				assert(hit_mol_sup.length() == 1);
				assert(hit_mol_sup.atEnd());
				const unique_ptr<ROMol> hit_mol_ptr(hit_mol_sup.next());
				auto& hitMol = *hit_mol_ptr;

				// Calculate Morgan fingerprint.
				const unique_ptr<SparseIntVect<uint32_t>> hitFp(getFingerprint(hitMol, 2));

				// Calculate Tanimoto similarity.
				const auto ts = TanimotoSimilarity(*qryFp, *hitFp);

				// Remove hydrogens to calculate canonical SMILES and descriptors.
				const unique_ptr<ROMol> hitMolNoH_ptr(removeHs(hitMol, false, false, false)); // implicitOnly, updateExplicitCount, sanitize
				const auto& hitMolNoH = *hitMolNoH_ptr;

				// Calculate canonical SMILES, molecular formula and descriptors. This can be done either by calculating them on the fly using the molecule with hydrogens removed, or by reading the precalculated values from *.u16 and *.f32 files.
				hitMol.setProp<unsigned int>("query", query_number);
				hitMol.setProp<double>("usrScore", usr1 ? u0score : u1score);
				hitMol.setProp<double>("usrcatScore", usr1 ? u1score : u0score);
				hitMol.setProp<double>("tanimotoScore", ts);
				hitMol.setProp<string>("database", cpdb.name);
				hitMol.setProp<string>("canonicalSMILES", MolToSmiles(hitMolNoH)); // Default parameters are: const ROMol& mol, bool doIsomericSmiles = true, bool doKekule = false, int rootedAtAtom = -1, bool canonical = true, bool allBondsExplicit = false, bool allHsExplicit = false, bool doRandom = false. https://www.rdkit.org/docs/cppapi/namespaceRDKit.html#a3636828cca83a233d7816f3652a9eb6b
				hitMol.setProp<string>("molFormula", calcMolFormula(hitMol)); // Calculate hydrogens in the molecular formula.
				hitMol.setProp<unsigned int>("numAtoms", cpdb.natm[k]);
				hitMol.setProp<unsigned int>("numHBD", cpdb.nhbd[k]);
				hitMol.setProp<unsigned int>("numHBA", cpdb.nhba[k]);
				hitMol.setProp<unsigned int>("numRotatableBonds", cpdb.nrtb[k]); // cpdb.nrtb[k] was precalculated from SMILES before adding hydrogens. Adding hydrogens may lead to more rotatable bonds. As a result, cpdb.nrtb[k] == calcNumRotatableBonds(hitMolNoH) <= calcNumRotatableBonds(hitMol)
				hitMol.setProp<unsigned int>("numRings", cpdb.nrng[k]);
				hitMol.setProp<double>("exactMW", cpdb.xmwt[k]);
				hitMol.setProp<double>("tPSA", cpdb.tpsa[k]);
				hitMol.setProp<double>("clogP", cpdb.clgp[k]);

				// Find heavy atoms.
				vector<vector<pair<int, int>>> matchVect;
				SubstructMatch(hitMol, *SubsetMols[0], matchVect);
				const auto num_matches = matchVect.size();
				assert(num_matches == hitMol.getNumHeavyAtoms());
				vector<int> hitHeavyAtoms(num_matches);
				for (size_t i = 0; i < num_matches; ++i)
				{
					hitHeavyAtoms[i] = matchVect[i].front().second;
//					assert(hitHeavyAtoms[i] == i); // Comment this assertion to avoid warning: comparison of integer expressions of different signedness. hitHeavyAtoms can be constructed using iota(hitHeavyAtoms.begin(), hitHeavyAtoms.end(), 0); because for RDKit-generated SDF compounds, heavy atom are always the first few atoms.
				}

				// Calculate the four reference points.
				const auto hitRefPoints = calcRefPoints(hitMol, hitHeavyAtoms);
				const Point3DConstPtrVect hitRefPointv
				{{
					&hitRefPoints[0],
					&hitRefPoints[1],
					&hitRefPoints[2],
					&hitRefPoints[3],
				}};

				// Calculate a 3D transform from the four reference points of the hit conformer to those of the query compound.
				Transform3D trans;
				AlignPoints(qryRefPointv, hitRefPointv, trans);

				// Apply the 3D transform to all atoms of the hit conformer.
				auto& hitCnf = hitMol.getConformer();
				transformConformer(hitCnf, trans);

				// Write the aligned hit conformer.
				hit_mol_sdf_writer.write(hitMol);
			}
			cout << local_time_string() << "Wrote " << hit_mol_sdf_per_qry_oss.tellp() << " bytes of hit molecules to a string stream" << endl;

			// If the size of the hitMolSdf field will not exceed 15MB after appending, then allow the appending. Reserve 1MB for the other fields, e.g. qryMolSdf.
			if (hit_mol_sdf_oss.tellp() + hit_mol_sdf_per_qry_oss.tellp() < 15000000)
			{
				hit_mol_sdf_oss << hit_mol_sdf_per_qry_oss.str();
				cout << local_time_string() << "Accumulated " << hit_mol_sdf_oss.tellp() << " bytes of hit molecules for the first " << ++query_number << " query molecules" << endl;
			}
			else
			{
				cout << local_time_string() << "Unable to accumulate " << hit_mol_sdf_per_qry_oss.tellp() << " bytes to the existing " << hit_mol_sdf_oss.tellp() << " bytes" << endl;
				break;
			}
		}
		const auto num_qry_mols_processed = query_number; // num_qry_mols_processed is the number of query molecules processed by the daemon.
		assert(num_qry_mols_processed <= num_qry_mols);
		cout << local_time_string() << "Processed " << num_qry_mols_processed << " out of " << num_qry_mols << " query molecules" << endl;

		// Update job status.
		const auto hit_mol_sdf = hit_mol_sdf_oss.str();
		cout << local_time_string() << "Writing " << hit_mol_sdf.size() << " bytes of hit molecules and setting end date" << endl;
		const auto endDate = system_clock::now();
		const int32_t cpdb_num_compounds = cpdb.num_compounds; // Create an int32_t intance, to be passed to kvp(). Caution: static_cast<int32_t>(cpdb.num_compounds) would cause the program to exit.
		const int32_t cpdb_num_conformers = cpdb.num_conformers;
		bsoncxx::builder::basic::document compt_update_builder;
		compt_update_builder.append(
			kvp("$set", [=](bsoncxx::builder::basic::sub_document set_subdoc) {
				set_subdoc.append(kvp("hitMolSdf", hit_mol_sdf));
				set_subdoc.append(kvp("endDate", bsoncxx::types::b_date(endDate)));
				set_subdoc.append(kvp("numQryMol", static_cast<int32_t>(num_qry_mols_processed))); // numQryMol is the number of query molecules completed by the daemon.
				set_subdoc.append(kvp("numLibMol", cpdb_num_compounds));
				set_subdoc.append(kvp("numLibCnf", cpdb_num_conformers));
				})
		);
		// http://mongocxx.org/api/current/classmongocxx_1_1collection.html#aece5216e5ae6fc3316c9da604f3b28f9
		const auto compt_update = coll.update_one(bsoncxx::builder::basic::make_document(kvp("_id", _id)), compt_update_builder.extract(), options::update()); // stdx::optional<result::update>. options: write_concern
		assert(compt_update);
		assert(compt_update->matched_count() == 1);
		assert(compt_update->modified_count() == 1);

		// Calculate runtime in seconds and screening speed in conformers per second.
		const auto runtime = duration_cast<nanoseconds>(endDate - startDate).count() * 1e-9; // Convert nanoseconds to seconds.
		const auto speed = cpdb.num_conformers * num_qry_mols_processed / runtime;
		cout
			<< local_time_string() << "Completed " << num_qry_mols_processed << " " << (num_qry_mols_processed == 1 ? "query" : "queries") << " in " << setprecision(3) << runtime << " seconds" << endl
			<< local_time_string() << "Screening speed was " << setprecision(0) << speed << " conformers per second" << endl
		;
	}
}
