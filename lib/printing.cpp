#pragma once
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>

namespace std {
template <>
struct is_integral<__int128> : std::integral_constant<bool, true> {};
} // namespace std

template <typename T>
std::ostream& operator<<(std::ostream& o, const std::optional<T>& v) {
	if(v) return o << *v;
	return o << "{}";
}

std::ostream& operator<<(std::ostream& o, const __int128& x) {
	if(x >= 1e18)
		return o << (long long)(x / (long long)1e18) << std::setw(18)
		         << (long long)((x < 0 ? -x : x) % (long long)1e18);
	else
		return o << (long long)x;
}
std::istream& operator>>(std::istream& o, __int128& x) {
	std::string input;
	o >> input;
	if(input.size() > 18) {
		std::string high = input.substr(0, input.size() - 18);
		std::string low  = input.substr(input.size() - 18, input.size());
		__int128 highval = stoull(high);
		__int128 lowval  = stoull(low);
		x                = highval * (long long)1e18 + (highval < 0 ? -1 : 1) * lowval;
	} else {
		x = stoull(input);
	}
	return o;
}
std::ostream& operator<<(std::ostream& o, const unsigned __int128& x) {
	if(x >= 1e18)
		return o << (long long)(x / (long long)1e18) << std::setw(18)
		         << (long long)(x % (long long)1e18);
	else
		return o << (long long)x;
}
std::istream& operator>>(std::istream& o, unsigned __int128& x) {
	std::string input;
	o >> input;
	if(input.size() > 18) {
		std::string high = input.substr(0, input.size() - 18);
		std::string low  = input.substr(input.size() - 18, input.size());
		x = __int128(stoull(high)) * __int128((long long)1e18) + __int128(stoull(low));
	} else {
		x = stoull(input);
	}
	return o;
}

template <typename X, typename Y>
std::ostream& operator<<(std::ostream& o, const std::pair<X, Y>& p) {
	return o << '(' << p.first << ',' << p.second << ')';
}
template <typename X, typename Y, typename Z>
std::ostream& operator<<(std::ostream& o, const std::tuple<X, Y, Z>& p) {
	return o << '(' << get<0>(p) << ',' << get<1>(p) << ',' << get<2>(p) << ')';
}
template <typename X, typename Y>
std::istream& operator>>(std::istream& i, std::pair<X, Y>& p) {
	return i >> p.first >> p.second;
}

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& o, const std::array<T, N>& v) {
	for(auto& x : v) o << x << "\t";
	return o;
}
template <typename T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
	for(const auto& x : v) o << x << "\t";
	return o;
}
template <typename U, typename V>
std::ostream& operator<<(std::ostream& o, const std::map<U, V>& m) {
	for(const auto& [k, v] : m) o << k << "\t -> " << v << "\n";
	return o;
}
template <typename T>
std::ostream& operator<<(std::ostream& o, const std::set<T>& m) {
	for(const auto& k : m) o << k << "\n";
	return o;
}
