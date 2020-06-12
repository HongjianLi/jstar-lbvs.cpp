#pragma once
#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <boost/date_time/posix_time/posix_time.hpp>
using namespace boost::posix_time;

inline static auto local_time()
{
	return to_simple_string(microsec_clock::local_time()) + " "; // TODO: update to std::chrono::format(std::chrono::system_clock::now()) when this c++20 feature is implemented in gcc or clang
}

#endif
