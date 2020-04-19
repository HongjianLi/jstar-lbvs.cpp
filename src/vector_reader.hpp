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
inline vector<T> read(const path src) // Sequential read can be very fast when using SSD
{
	ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	cout << local_time() << "Reading " << src.filename() << " of " << num_bytes << " bytes" << endl;
	vector<T> buf;
	buf.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(buf.data()), num_bytes);
	return buf;
}

inline vector<string> readLines(const path src) // Sequential read can be very fast when using SSD
{
	ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	cout << local_time() << "Reading " << src.filename() << " of " << num_bytes << " bytes" << endl;
	vector<string> vec;
	vec.reserve(num_bytes / 6); // An average line size is 6 bytes.
	ifs.seekg(0);
	string line;
	while (getline(ifs, line)) {
		vec.push_back(move(line));
	}
	return vec;
}

//! Represents a vector of headers, which are the positions of the first character of each entity in a vector file.
template <typename size_type> // Size type of the encoded binary representation of the footer file. This is usually size_t in order to support seeking vector files of >4GB.
class header_vector
{
public:
	//! Reads the footer file for the vector file.
	explicit header_vector(path src) :
		ftr(read<size_t>(src.replace_extension(src.extension().string() + ".ftr")))
	{
	}

	size_t size() const
	{
		return ftr.size();
	}

protected:
	vector<size_type> ftr;
};

template <typename size_type>
class stream_vector : public header_vector<size_type>
{
public:
	explicit stream_vector(const path src) :
		header_vector<size_type>(src),
		ifs(src, ios::binary)
	{
	}

	string operator[](const size_t index)
	{
		const auto pos = index ? this->ftr[index - 1] : 0;
		const auto len = this->ftr[index] - pos;
		string buf;
		buf.resize(len);
		ifs.seekg(pos);
		ifs.read(const_cast<char*>(buf.data()), len);
		return buf;
	}

protected:
	ifstream ifs;
};

#endif
