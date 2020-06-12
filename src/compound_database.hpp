#pragma once
#ifndef COMPOUND_DATABASE_HPP
#define COMPOUND_DATABASE_HPP

#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
using namespace std;
using namespace std::filesystem;

//! Represents a compound database.
class compound_database
{
public:
	explicit compound_database(const path dpth);
	string read_conformer(const size_t index, ifstream& ifs) const;

	string name; //!< Database name.
	path dpth; //!< Path to the database directory.
	size_t num_compounds; //!< Number of compound.
	size_t num_conformers; //!< Number of 3D conformers. Usually num_conformers == num_compounds < 2
	vector<string> cpid; //!< Compound identifies.
	vector<uint16_t> natm; //!< Number of atoms.
	vector<uint16_t> nhbd; //!< Number of hydrogen bond donors.
	vector<uint16_t> nhba; //!< Number of hydrogen bond acceptors.
	vector<uint16_t> nrtb; //!< Number of rotatable bonds.
	vector<uint16_t> nrng; //!< Number of rings.
	vector<float> xmwt; //!< Exact molecular weight.
	vector<float> tpsa; //!< Topological polar surface area.
	vector<float> clgp; //!< clogP
	vector<array<float, 60>> usrcat; //!< USRCAT features.
	vector<size_t> conformers_sdf_ftr; //!< Footer file of conformers.sdf
protected:
	template <typename T>
	static void read_types(const path src, vector<T>& vec); // Sequential read can be very fast when using SSD.
	static void read_lines(const path src, vector<string>& vec); // Sequential read can be very fast when using SSD.
	static string read_string(const vector<size_t>& ftr, const size_t index, ifstream& ifs);
};

#endif
