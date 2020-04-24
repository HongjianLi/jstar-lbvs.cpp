#pragma once
#ifndef VECTOR_READER_HPP
#define VECTOR_READER_HPP

#include <vector>
#include <filesystem>
#include <boost/date_time/posix_time/posix_time.hpp>
using namespace std;
using namespace std::filesystem;
using namespace boost::posix_time;

inline static auto local_time()
{
	return to_simple_string(microsec_clock::local_time()) + " "; // TODO: update to std::chrono::format(std::chrono::system_clock::now()) when this c++20 feature is implemented in gcc or clang
}

template <typename T>
inline void read_types(const path src, vector<T>& vec) // Sequential read can be very fast when using SSD.
{
	ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	cout << local_time() << "Reading " << src.filename() << " of " << num_bytes << " bytes" << endl;
	vec.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(vec.data()), num_bytes);
}

inline void read_lines(const path src, vector<string>& vec) // Sequential read can be very fast when using SSD.
{
	ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	cout << local_time() << "Reading " << src.filename() << " of " << num_bytes << " bytes" << endl;
	vec.reserve(num_bytes / 17); // Assuming an average line size is 17 bytes (e.g. a ZINC ID line has 16 characters plus \n).
	ifs.seekg(0);
	string line;
	while (getline(ifs, line)) {
		vec.push_back(move(line));
	}
}

inline string read_string(const vector<size_t>& ftr, const size_t index, ifstream& ifs)
{
	const auto pos = index ? ftr[index - 1] : 0;
	const auto len = ftr[index] - pos;
	string str;
	str.resize(len);
	ifs.seekg(pos);
	ifs.read(const_cast<char*>(str.data()), len);
	return str;
}

#endif
