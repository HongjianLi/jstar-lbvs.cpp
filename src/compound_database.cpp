#include <iostream>
#include "utility.hpp"
#include "compound_database.hpp"

compound_database::compound_database(const path dpth) : dpth(dpth)
{
	// Obtain database name.
	name = dpth.filename().string();
	cout << local_time_string() << "Reading " << name << endl;

	// Read molecular descriptor files.
	read_types<uint16_t>(dpth / "natm.u16", natm);
	num_compounds = natm.size();
	read_types<uint16_t>(dpth / "nhbd.u16", nhbd);
	assert(nhbd.size() == num_compounds);
	read_types<uint16_t>(dpth / "nhba.u16", nhba);
	assert(nhba.size() == num_compounds);
	read_types<uint16_t>(dpth / "nrtb.u16", nrtb);
	assert(nrtb.size() == num_compounds);
	read_types<uint16_t>(dpth / "nrng.u16", nrng);
	assert(nrng.size() == num_compounds);
	read_types<float>(dpth / "xmwt.f32", xmwt);
	assert(xmwt.size() == num_compounds);
	read_types<float>(dpth / "tpsa.f32", tpsa);
	assert(tpsa.size() == num_compounds);
	read_types<float>(dpth / "clgp.f32", clgp);
	assert(clgp.size() == num_compounds);

	// Read usrcat feature file.
	read_types<array<float, 60>>(dpth / "usrcat.f32", usrcat);
	num_conformers = usrcat.size();
	assert(num_conformers == num_compounds << 2);

	// Read conformers.sdf footer file.
	read_types<size_t>(dpth / "conformers.sdf.ftr.u64", conformers_sdf_ftr);
	assert(conformers_sdf_ftr.size() == num_conformers);

	cout << local_time_string() << "Found " << num_compounds << " compounds and " << num_conformers << " conformers" << endl;
};

string compound_database::read_conformer(const size_t index, ifstream& ifs) const
{
	return read_string(conformers_sdf_ftr, index, ifs);
}

template <typename T>
void compound_database::read_types(const path src, vector<T>& vec) // Sequential read can be very fast when using SSD.
{
	ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	cout << local_time_string() << "Reading " << src.filename() << " of " << num_bytes << " bytes" << endl;
	vec.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(vec.data()), num_bytes);
}

// Read a string from the ifstream after seeking to the corresponding position obtained from the footer vector.
string compound_database::read_string(const vector<size_t>& ftr, const size_t index, ifstream& ifs)
{
	const auto pos = index ? ftr[index - 1] : 0;
	const auto len = ftr[index] - pos;
	string str;
	str.resize(len);
	ifs.seekg(pos);
	ifs.read(str.data(), len);
	return str;
}
