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
#include <boost/algorithm/string.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <mongocxx/instance.hpp>
#include <mongocxx/pool.hpp>
#include <mongocxx/client.hpp>
#include <bsoncxx/json.hpp>
#include "io_service_pool.hpp"
#include "safe_counter.hpp"
#include "compound_database.hpp"
#include "vector_reader.hpp"
using namespace std;
using namespace std::chrono;
using namespace std::filesystem;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDDepict;
using namespace RDKit::MorganFingerprints;
using namespace RDGeom;
using namespace RDNumeric::Alignments;
using namespace MolTransforms;
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
		cout << "lbvs host port user pwd dbs_path" << endl;
		return 0;
	}

	// Fetch command line arguments.
	const auto host = argv[1];
	const auto port = argv[2];
	const auto user = argv[3];
	const auto pwd = argv[4];
	const path cpdbs_path = argv[5];

	// Initialize constants.
	cout << local_time() << "Initializing" << endl;
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
	cout << cpdbs_path << endl;
	for (directory_iterator cpdbs_dir_iter(cpdbs_path), end_dir_iter; cpdbs_dir_iter != end_dir_iter; ++cpdbs_dir_iter)
	{
		// Filter out non directory.
		if (!cpdbs_dir_iter->is_directory()) continue;

		// Create a compound database instance.
		databases.emplace_back();
		auto& cpdb = databases.back();

		// Assign database path and name.
		cpdb.dpth = cpdbs_dir_iter->path();
		cpdb.name = cpdb.dpth.filename().string();

		// Read id file.
		cout << local_time() << "Reading " << cpdb.name << endl;
		read_lines(cpdb.dpth / "id.txt", cpdb.cpid);
		cpdb.num_compounds = cpdb.cpid.size();
		cout << local_time() << "Found " << cpdb.num_compounds << " compounds" << endl;

		// Read molecular descriptor files.
		read_types<uint16_t>(cpdb.dpth / "natm.u16", cpdb.natm);
		assert(cpdb.natm.size() == cpdb.num_compounds);
		read_types<uint16_t>(cpdb.dpth / "nhbd.u16", cpdb.nhbd);
		assert(cpdb.nhbd.size() == cpdb.num_compounds);
		read_types<uint16_t>(cpdb.dpth / "nhba.u16", cpdb.nhba);
		assert(cpdb.nhba.size() == cpdb.num_compounds);
		read_types<uint16_t>(cpdb.dpth / "nrtb.u16", cpdb.nrtb);
		assert(cpdb.nrtb.size() == cpdb.num_compounds);
		read_types<uint16_t>(cpdb.dpth / "nrng.u16", cpdb.nrng);
		assert(cpdb.nrng.size() == cpdb.num_compounds);
		read_types<float>(cpdb.dpth / "xmwt.f32", cpdb.xmwt);
		assert(cpdb.xmwt.size() == cpdb.num_compounds);
		read_types<float>(cpdb.dpth / "tpsa.f32", cpdb.tpsa);
		assert(cpdb.tpsa.size() == cpdb.num_compounds);
		read_types<float>(cpdb.dpth / "clgp.f32", cpdb.clgp);
		assert(cpdb.clgp.size() == cpdb.num_compounds);

		// Read usrcat feature file.
		read_types<array<float, 60>>(cpdb.dpth / "usrcat.f32", cpdb.usrcat);
		cpdb.num_conformers = cpdb.usrcat.size();
		cout << local_time() << "Found " << cpdb.num_conformers << " conformers" << endl;
		assert(cpdb.num_conformers == cpdb.num_compounds << 2);

		// Read conformers.sdf and descriptors.tsv footer files.
		read_types<size_t>(cpdb.dpth / "conformers.sdf.ftr", cpdb.conformers_sdf_ftr);
		assert(cpdb.conformers_sdf_ftr.size() == cpdb.num_conformers);
		read_types<size_t>(cpdb.dpth / "descriptors.tsv.ftr", cpdb.descriptors_tsv_ftr);
		assert(cpdb.descriptors_tsv_ftr.size() == cpdb.num_compounds);
	}
	string db_op_in_array;
	for (size_t i = 0; i < databases.size(); ++i)
	{
		db_op_in_array += (i ? ",\"" : "\"") + databases[i].name + "\"";
	}

	// Connect to mongodb and authenticate user.
	cout << local_time() << "Connecting to " << host << ':' << port << " and authenticating " << user << endl;
	const instance inst; // The constructor and destructor initialize and shut down the driver. http://mongocxx.org/api/current/classmongocxx_1_1instance.html
	const uri uri("mongodb://localhost:27017/?minPoolSize=0&maxPoolSize=2"); // When connecting to a replica set, it is much more efficient to use a pool as opposed to manually constructing client objects.
	pool pool(uri);
	const auto client = pool.acquire(); // Return value of acquire() is an instance of entry. An entry is a handle on a client object acquired via the pool.
	auto coll = client->database("jstar").collection("lbvs");
	const auto jobid_filter = bsoncxx::from_json(R"({ "startDate" : { "$exists" : false }, "database": { "$in": [)" + db_op_in_array + R"(] }})");
	const auto jobid_foau_options = options::find_one_and_update().sort(bsoncxx::from_json(R"({ "submitDate" : 1 })")).projection(bsoncxx::from_json(R"({ "_id" : 1, "qryMolSdf": 1, "database": 1, "score": 1 })")); // By default, the original document is returned

	// Initialize variables.
	array<vector<int>, num_subsets> subsets;
	array<vector<double>, num_refPoints> dista;
	alignas(32) array<double, 60> q;

	// Initialize an io service pool and create worker threads for later use.
	const size_t num_threads = thread::hardware_concurrency();
	cout << local_time() << "Creating an io service pool of " << num_threads << " worker threads" << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;
	const auto num_chunks = num_threads << 2;

	// Enter event loop.
	cout << local_time() << "Entering event loop" << endl;
	cout.setf(ios::fixed, ios::floatfield);
	bool sleeping = false;
	while (true)
	{
		// Fetch an incompleted job in a first-come-first-served manner.
		if (!sleeping) cout << local_time() << "Fetching an incompleted job" << endl;
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
			if (!sleeping) cout << local_time() << "Sleeping" << endl;
			sleeping = true;
			this_thread::sleep_for(std::chrono::seconds(2));
			continue;
		}
		sleeping = false;
		const auto jobid_view = jobid_document->view();

		// Obtain job properties.
		const auto _id = jobid_view["_id"].get_oid().value;
		cout << local_time() << "Executing job " << _id.to_string() << endl;
		const auto qry_mol_sdf = jobid_view["qryMolSdf"].get_utf8().value; // get_utf8().value returns an instance of std::string_view.
		const auto cpdb_name = jobid_view["database"].get_utf8().value;
		const auto score = jobid_view["score"].get_utf8().value;
		assert(score.compare("USR") == 0 || score.compare("USRCAT") == 0);
		const size_t usr0 = score.compare("USRCAT") == 0; // Specify the primary sorting score. 0: USR; 1: USRCAT.
		const auto usr1 = usr0 ^ 1;
		const auto qnu0 = qn[usr0];
		const auto qnu1 = qn[usr1];

		// Obtain a constant reference to the selected database.
		cout << local_time() << "Finding the selected compound database" << endl;
		size_t cpdb_index = 0;
		while (cpdb_name.compare(databases[cpdb_index].name) != 0) ++cpdb_index;
		const auto& cpdb = databases[cpdb_index];

		// Read the user-supplied SDF file.
		cout << local_time() << "Reading the query file" << endl;
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
		const auto chunk_size = 1 + (cpdb.num_compounds - 1) / num_chunks;
		assert(chunk_size * num_chunks >= cpdb.num_compounds);
		assert(chunk_size >= num_hits);
		cout << local_time() << "Using " << num_chunks << " chunks and a chunk size of " << chunk_size << endl;
		vector<size_t> scase(cpdb.num_compounds);
		vector<size_t> zcase(num_hits * (num_chunks - 1) + min(num_hits, cpdb.num_compounds - chunk_size * (num_chunks - 1))); // The last chunk might have fewer than num_hits records.
		ifstream conformers_sdf_ifs(cpdb.dpth / "conformers.sdf");
		ifstream descriptors_tsv_ifs(cpdb.dpth / "descriptors.tsv");

		// Process each of the query compounds sequentially.
		ostringstream hit_mol_sdf_oss, hit_mol_csv_oss;
		SDWriter hit_mol_sdf_writer(&hit_mol_sdf_oss, false); // std::ostream*, bool takeOwnership
		const auto num_queries = 1; // Restrict the number of query compounds to 1. Setting num_queries = qry_mol_sup.length() to execute any number of query compounds.
		for (unsigned int query_number = 0; query_number < num_queries; ++query_number)
		{
			cout << local_time() << "Parsing query compound " << query_number << endl;
			const unique_ptr<ROMol> qry_mol_ptr(qry_mol_sup.next()); // Calling next() may print "ERROR: Could not sanitize compound on line XXXX" to stderr.
			auto& qryMol = *qry_mol_ptr;

			// Get the number of atoms, including and excluding hydrogens.
			const auto num_atoms = qryMol.getNumAtoms();
			const auto num_heavy_atoms = qryMol.getNumHeavyAtoms();
			assert(num_heavy_atoms);
			cout << local_time() << "Found " << num_atoms << " atoms and " << num_heavy_atoms << " heavy atoms" << endl;

			// Calculate Morgan fingerprint.
			cout << local_time() << "Calculating Morgan fingerprint" << endl;
			const unique_ptr<SparseIntVect<uint32_t>> qryFp(getFingerprint(qryMol, 2));

			// Classify atoms to pharmacophoric subsets.
			cout << local_time() << "Classifying atoms into subsets" << endl;
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
				cout << local_time() << "Found " << num_matches << " atoms for subset " << k << endl;
			}
			const auto& subset0 = subsets.front();
			assert(subset0.size() == num_heavy_atoms);

			// Calculate the four reference points.
			cout << local_time() << "Calculating " << num_refPoints << " reference points" << endl;
			const auto qryRefPoints = calcRefPoints(qryMol, subset0);
			const Point3DConstPtrVect qryRefPointv
			{{
				&qryRefPoints[0],
				&qryRefPoints[1],
				&qryRefPoints[2],
				&qryRefPoints[3],
			}};

			// Precalculate the distances of heavy atoms to the reference points, given that subsets[1 to 4] are subsets of subsets[0].
			cout << local_time() << "Calculating " << num_heavy_atoms * num_refPoints << " pairwise distances" << endl;
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
			cout << local_time() << "Calculating " << 3 * num_refPoints * num_subsets << " moments of USRCAT feature" << endl;
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
			cout << local_time() << "Screening " << cpdb.name << " and calculating " << cpdb.num_compounds << " " << usr_names[usr0] << " scores from " << cpdb.num_conformers << " conformers" << endl;
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
			cout << local_time() << "Sorting " << zcase.size() << " hits by " << usr_names[usr0] << " score" << endl;
			sort(zcase.begin(), zcase.end(), compare);

			// Create output directory and write output files.
			cout << local_time() << "Writing output string streams" << endl;
			hit_mol_csv_oss.setf(ios::fixed, ios::floatfield);
			hit_mol_csv_oss << setprecision(8) << "ID,USR score,USRCAT score,2D Tanimoto score,canonicalSMILES,molFormula,natm,nhbd,nhba,nrtb,nrng,xmwt,tpsa,clgp\n";
			for (size_t l = 0; l < num_hits; ++l)
			{
				// Obtain indexes to the hit compound and the hit conformer.
				const auto k = zcase[l];
				const auto j = cnfids[k];

				// Read SDF content of the hit conformer.
				istringstream hit_mol_sdf_iss(read_string(cpdb.conformers_sdf_ftr, j, conformers_sdf_ifs));

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

				// Find heavy atoms.
				vector<vector<pair<int, int>>> matchVect;
				SubstructMatch(hitMol, *SubsetMols[0], matchVect);
				const auto num_matches = matchVect.size();
				assert(num_matches == hitMol.getNumHeavyAtoms());
				vector<int> hitHeavyAtoms(num_matches);
				for (size_t i = 0; i < num_matches; ++i)
				{
					hitHeavyAtoms[i] = matchVect[i].front().second;
					assert(hitHeavyAtoms[i] == i); // hitHeavyAtoms can be constructed using iota(hitHeavyAtoms.begin(), hitHeavyAtoms.end(), 0); because for RDKit-generated SDF compounds, heavy atom are always the first few atoms.
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

				// Calculate the secondary score of the saved conformer, which has the best primary score.
				const auto& d = cpdb.usrcat[j];
				double s = 0;
				for (size_t i = 0; i < qnu1; ++i)
				{
					s += abs(q[i] - d[i]);
				}
				const auto u0score = 1 / (1 + scores[k] * qv[usr0]); // Primary score of the current compound.
				const auto u1score = 1 / (1 + s         * qv[usr1]); // Secondary score of the current compound.

				// Read SDF content of the hit conformer.
				vector<string> descriptors;
				split(descriptors, read_string(cpdb.descriptors_tsv_ftr, k, descriptors_tsv_ifs), boost::is_any_of("	")); // Split the descriptor line into columns, which are [ID	canonicalSMILES	molFormula	natm	nhbd	nhba	nrtb	nrng	xmwt	tpsa	clgp]
				assert(descriptors[0] == cpdb.cpid[k]);
				assert(stoul(descriptors[3]) == cpdb.natm[k]);
				assert(stoul(descriptors[4]) == cpdb.nhbd[k]);
				assert(stoul(descriptors[5]) == cpdb.nhba[k]);
				assert(stoul(descriptors[6]) == cpdb.nrtb[k]);
				assert(stoul(descriptors[7]) == cpdb.nrng[k]);
				assert(stof(descriptors[8]) == cpdb.xmwt[k]);
				assert(stof(descriptors[9]) == cpdb.tpsa[k]);
				assert(stof(descriptors[10]) == cpdb.clgp[k]);
				hit_mol_csv_oss
					<< cpdb.cpid[k] // ID
//					<< ',' << cpdb.name
					<< ',' << (usr1 ? u0score : u1score)
					<< ',' << (usr1 ? u1score : u0score)
					<< ',' << ts // Tanimoto score
					<< ',' << descriptors[1] // Canonical SMILES
					<< ',' << descriptors[2] // Molecular formula
					<< ',' << cpdb.natm[k]
					<< ',' << cpdb.nhbd[k]
					<< ',' << cpdb.nhba[k]
					<< ',' << cpdb.nrtb[k]
					<< ',' << cpdb.nrng[k]
					<< ',' << cpdb.xmwt[k]
					<< ',' << cpdb.tpsa[k]
					<< ',' << cpdb.clgp[k]
					<< '\n'
				;
			}
		}

		// Update job status.
		cout << local_time() << "Setting end date" << endl;
		const auto endDate = system_clock::now();
		const auto hit_mol_sdf = hit_mol_sdf_oss.str();
		const auto hit_mol_csv = hit_mol_csv_oss.str();
		bsoncxx::builder::basic::document compt_update_builder;
		compt_update_builder.append(
			kvp("$set", [=](bsoncxx::builder::basic::sub_document set_subdoc) {
				set_subdoc.append(kvp("hitMolSdf", hit_mol_sdf));
				set_subdoc.append(kvp("hitMolCsv", hit_mol_csv));
				set_subdoc.append(kvp("endDate", bsoncxx::types::b_date(endDate)));
				set_subdoc.append(kvp("numQueries", num_queries));
				set_subdoc.append(kvp("numConformers", static_cast<int64_t>(cpdb.num_conformers)));
			})
		);
		// http://mongocxx.org/api/current/classmongocxx_1_1collection.html#aece5216e5ae6fc3316c9da604f3b28f9
		const auto compt_update = coll.update_one(bsoncxx::builder::basic::make_document(kvp("_id", _id)), compt_update_builder.extract(), options::update()); // stdx::optional<result::update>. options: write_concern
		assert(compt_update);
		assert(compt_update->matched_count() == 1);
		assert(compt_update->modified_count() == 1);

		// Calculate runtime in seconds and screening speed in thousand conformers per second.
		const auto runtime = (endDate - startDate).count() * 1e-9; // in seconds
		const auto speed = cpdb.num_conformers * num_queries * 1e-3 / runtime;
		cout
			<< local_time() << "Completed " << num_queries << " " << (num_queries == 1 ? "query" : "queries") << " in " << setprecision(3) << runtime << " seconds" << endl
			<< local_time() << "Screening speed was " << setprecision(2) << speed << " K conformers per second" << endl
		;
	}
}
