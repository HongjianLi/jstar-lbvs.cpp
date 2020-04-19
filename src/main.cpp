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
#include "vector_reader.hpp"
#include "safe_counter.hpp"
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
	ctd /= num_points;
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
		cout << "lbvs host port user pwd jobs_path" << endl;
		return 0;
	}

	// Fetch command line arguments.
	const auto host = argv[1];
	const auto port = argv[2];
	const auto user = argv[3];
	const auto pwd = argv[4];
	const path jobs_path = argv[5];

	// Connect to host and authenticate user.
	cout << local_time() << "Connecting to " << host << ':' << port << " and authenticating " << user << endl;
	const instance inst; // The constructor and destructor initialize and shut down the driver. http://mongocxx.org/api/current/classmongocxx_1_1instance.html
	const uri uri("mongodb://localhost:27017/?minPoolSize=0&maxPoolSize=2"); // When connecting to a replica set, it is much more efficient to use a pool as opposed to manually constructing client objects.
	pool pool(uri);
	const auto client = pool.acquire(); // Return value of acquire() is an instance of entry. An entry is a handle on a client object acquired via the pool.
	const auto db = client->database("jstar");
	auto coll = db.collection("usr2");
	const auto jobid_filter = bsoncxx::from_json(R"({ "started" : { "$exists" : false }})");
	const auto jobid_foau_options = options::find_one_and_update().sort(bsoncxx::from_json(R"({ "submitted" : 1 })")).projection(bsoncxx::from_json(R"({ "_id" : 1, "usr": 1 })")); // By default, the original document is returned

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

	// Read id file.
	const path dbPath = "databases";
	const string collName = "Selleckchem";
	const path collPath = dbPath / collName;
	cout << local_time() << "Reading " << collName << endl;
//	const auto id_u32 = read<uint32_t>(collPath / "id.u32");
	const auto id_str = readLines(collPath / "id.txt");
//	const auto num_compounds = id_u32.size();
	const auto num_compounds = id_str.size();
	cout << local_time() << "Found " << num_compounds << " compounds from " << collName << endl;

/*	// Read property files.
	const auto numAtoms_u16 = read<uint16_t>(collPath / "numAtoms.u16");
	assert(numAtoms_u16.size() == num_compounds);
	const auto numHBD_u16 = read<uint16_t>(collPath / "numHBD.u16");
	assert(numHBD_u16.size() == num_compounds);
	const auto numHBA_u16 = read<uint16_t>(collPath / "numHBA.u16");
	assert(numHBA_u16.size() == num_compounds);
	const auto numRotatableBonds_u16 = read<uint16_t>(collPath / "numRotatableBonds.u16");
	assert(numRotatableBonds_u16.size() == num_compounds);
	const auto numRings_u16 = read<uint16_t>(collPath / "numRings.u16");
	assert(numRings_u16.size() == num_compounds);
	const auto exactMW_f32 = read<float>(collPath / "exactMW.f32");
	assert(exactMW_f32.size() == num_compounds);
	const auto tPSA_f32 = read<float>(collPath / "tPSA.f32");
	assert(tPSA_f32.size() == num_compounds);
	const auto clogP_f32 = read<float>(collPath / "clogP.f32");
	assert(clogP_f32.size() == num_compounds);*/

	// Read usrcat feature file.
	const auto usrcat_f32 = read<array<float, qn.back()>>(collPath / "usrcat.f32");
	const auto num_conformers = usrcat_f32.size();
	cout << local_time() << "Found " << num_conformers << " conformers from " << collName << endl;
	assert(num_conformers == num_compounds << 2);

	// Read ligand footer file and open ligand SDF file for seeking and reading.
	stream_vector<size_t> descriptors(collPath / "descriptors.tsv");
	assert(descriptors.size() == num_compounds);
	stream_vector<size_t> conformers(collPath / "conformers.sdf");
	assert(conformers.size() == num_conformers);

	// Initialize variables.
	array<vector<int>, num_subsets> subsets;
	array<vector<double>, num_refPoints> dista;
	alignas(32) array<double, qn.back()> q;

	// Initialize vectors to store compounds' primary score and their corresponding conformer.
	vector<double> scores(num_compounds); // Primary score of compounds.
	vector<size_t> cnfids(num_compounds); // ID of conformer with the best primary score.
	const auto compare = [&](const size_t val0, const size_t val1) // Sort by the primary score.
	{
		return scores[val0] < scores[val1];
	};

	// Initialize an io service pool and create worker threads for later use.
	const size_t num_threads = thread::hardware_concurrency();
	cout << local_time() << "Creating an io service pool of " << num_threads << " worker threads" << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;

	// Initialize the number of chunks and the number of compounds per chunk. TODO: the choice of num_chunks depends on num_compounds.
	const auto num_chunks = num_threads << 2;
	const auto chunk_size = 1 + (num_compounds - 1) / num_chunks;
	assert(chunk_size * num_chunks >= num_compounds);
	assert(chunk_size >= num_hits);
	cout << local_time() << "Using " << num_chunks << " chunks and a chunk size of " << chunk_size << endl;
	vector<size_t> scase(num_compounds);
	vector<size_t> zcase(num_hits * (num_chunks - 1) + min(num_hits, num_compounds - chunk_size * (num_chunks - 1))); // The last chunk might have fewer than num_hits records.

	// Enter event loop.
	cout << local_time() << "Entering event loop" << endl;
	cout.setf(ios::fixed, ios::floatfield);
	bool sleeping = false;
	while (true)
	{
		// Fetch an incompleted job in a first-come-first-served manner.
		if (!sleeping) cout << local_time() << "Fetching an incompleted job" << endl;
		const auto started = system_clock::now();
		bsoncxx::builder::basic::document jobid_update_builder;
		jobid_update_builder.append(
			kvp("$set", [=](bsoncxx::builder::basic::sub_document set_subdoc) {
				set_subdoc.append(kvp("started", bsoncxx::types::b_date(started)));
			})
		);
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
		const auto job_path = jobs_path / _id.to_string();
		const size_t usr0 = jobid_view["usr"].get_int32(); // Specify the primary sorting score. 0: USR; 1: USRCAT.
		assert(usr0 == 0 || usr0 == 1);
		const auto usr1 = usr0 ^ 1;
		const auto qnu0 = qn[usr0];
		const auto qnu1 = qn[usr1];

		// Read the user-supplied SDF file.
		cout << local_time() << "Reading the query file" << endl;
		SDMolSupplier sup((job_path / "query.sdf").string(), true, false, true); // sanitize, removeHs, strictParsing. Note: setting removeHs=true (which is the default setting) will lead to fewer hydrogen bond acceptors being matched.

		// Process each of the query compounds sequentially.
		const auto num_queries = 1; // Restrict the number of query compounds to 1. Setting num_queries = sup.length() to execute any number of query compounds.
		for (unsigned int query_number = 0; query_number < num_queries; ++query_number)
		{
			cout << local_time() << "Parsing query compound " << query_number << endl;
			const unique_ptr<ROMol> qry_ptr(sup.next()); // Calling next() may print "ERROR: Could not sanitize compound on line XXXX" to stderr.
			auto& qryMol = *qry_ptr;

			// Get the number of atoms, including and excluding hydrogens.
			const auto num_atoms = qryMol.getNumAtoms();
			const auto num_heavy_atoms = qryMol.getNumHeavyAtoms();
			assert(num_heavy_atoms);
			cout << local_time() << "Found " << num_atoms << " atoms and " << num_heavy_atoms << " heavy atoms" << endl;

			// Create an output directory.
			cout << local_time() << "Creating output directory" << endl;
			const auto output_dir = job_path / to_string(query_number);
			create_directory(output_dir);

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
			cout << local_time() << "Calculating " << num_compounds << " " << usr_names[usr0] << " scores" << endl;
			scores.assign(scores.size(), numeric_limits<double>::max());
			iota(scase.begin(), scase.end(), 0);
			cnt.init(num_chunks);
			for (size_t l = 0; l < num_chunks; ++l)
			{
				io.post([&,l]()
				{
					// Loop over compounds of the current chunk.
					const auto chunk_beg = chunk_size * l;
					const auto chunk_end = min(chunk_beg + chunk_size, num_compounds);
					for (size_t k = chunk_beg; k < chunk_end; ++k)
					{
						// Loop over conformers of the current compound and calculate their primary score.
						auto& scorek = scores[k];
						for (size_t j = k << 2; j < (k + 1) << 2; ++j)
						{
							const auto& d = usrcat_f32[j];
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
			cout << local_time() << "Writing output files" << endl;
			SDWriter hits_sdf((output_dir / "hits.sdf").string());
			ofstream hits_csv(output_dir / "hits.csv");
			hits_csv.setf(ios::fixed, ios::floatfield);
			hits_csv << setprecision(8) << "ID,USR score,USRCAT score,2D Tanimoto score,canonicalSMILES,molFormula,numAtoms,numHBD,numHBA,numRotatableBonds,numRings,exactMW,tPSA,clogP\n";
			for (size_t l = 0; l < num_hits; ++l)
			{
				// Obtain indexes to the hit compound and the hit conformer.
				const auto k = zcase[l];
				const auto j = cnfids[k];

				// Read SDF content of the hit conformer.
				const auto hitSdf = conformers[j];

				// Construct a RDKit ROMol object.
				istringstream iss(hitSdf);
				SDMolSupplier sup(&iss, false, true, false, true);
				assert(sup.length() == 1);
				assert(sup.atEnd());
				const unique_ptr<ROMol> hit_ptr(sup.next());
				auto& hitMol = *hit_ptr;

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
				hits_sdf.write(hitMol);

				// Calculate the secondary score of the saved conformer, which has the best primary score.
				const auto& d = usrcat_f32[j];
				double s = 0;
				for (size_t i = 0; i < qnu1; ++i)
				{
					s += abs(q[i] - d[i]);
				}

				const auto u0score = 1 / (1 + scores[k] * qv[usr0]); // Primary score of the current compound.
				const auto u1score = 1 / (1 + s         * qv[usr1]); // Secondary score of the current compound.
				const auto id = id_str[k];
				vector<string> descs;
				split(descs, descriptors[k], boost::is_any_of("	")); // Split the descriptor line into columns, which are [ID	canonicalSMILES	molFormula	numAtoms	numHBD	numHBA	numRotatableBonds	numRings	exactMW	tPSA	clogP	subset]
				hits_csv
//					<< boost::format("%08d") % id // SCUBIDOO
//					<< boost::format("S%d") % id // Selleckchem
					<< id
//					<< ',' << descs[1]
//					<< ',' << collName
					<< ',' << (usr1 ? u0score : u1score)
					<< ',' << (usr1 ? u1score : u0score)
					<< ',' << ts
//					<< ',' // subset
					<< ',' << descs[1]
					<< ',' << descs[2]
					<< ',' << descs[3]
					<< ',' << descs[4]
					<< ',' << descs[5]
					<< ',' << descs[6]
					<< ',' << descs[7]
					<< ',' << descs[8]
					<< ',' << descs[9]
					<< ',' << descs[10]
					<< '\n'
				;
			}
		}

		// Update job status.
		cout << local_time() << "Setting completed time" << endl;
		const auto completed = system_clock::now();
		bsoncxx::builder::basic::document compt_update_builder;
		compt_update_builder.append(
			kvp("$set", [=](bsoncxx::builder::basic::sub_document set_subdoc) {
				set_subdoc.append(kvp("completed", bsoncxx::types::b_date(completed)));
				set_subdoc.append(kvp("nqueries", num_queries));
				set_subdoc.append(kvp("numConformers", static_cast<int64_t>(num_conformers)));
			})
		);
		const auto compt_update = coll.update_one(bsoncxx::builder::basic::make_document(kvp("_id", _id)), compt_update_builder.extract(), options::update()); // stdx::optional<result::update>. options: write_concern
		assert(compt_update);
		assert(compt_update->matched_count() == 1);
		assert(compt_update->modified_count() == 1);

		// Calculate runtime in seconds and screening speed in million conformers per second.
		const auto runtime = (completed - started).count() * 1e-9; // in seconds
		const auto speed = num_conformers * 1e-6 * num_queries / runtime;
		cout
			<< local_time() << "Completed " << num_queries << " " << (num_queries == 1 ? "query" : "queries") << " in " << setprecision(3) << runtime << " seconds" << endl
			<< local_time() << "Screening speed was " << setprecision(0) << speed << " M conformers per second" << endl
		;
	}
}
