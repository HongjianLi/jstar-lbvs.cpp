#pragma once
#ifndef COMPOUND_DATABASE_HPP
#define COMPOUND_DATABASE_HPP

#include <string>
#include <vector>
#include <filesystem>
using namespace std;
using namespace std::filesystem;

//! Represents a compound database.
class compound_database
{
public:
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
	vector<size_t> descriptors_tsv_ftr; //!< Footer file of descriptors.tsv
};

#endif
