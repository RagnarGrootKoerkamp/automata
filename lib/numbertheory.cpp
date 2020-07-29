#pragma once

#include <algorithm>
#include <cassert>

template <typename T>
T gcd(T a, T b) {
	while(b != T{0}) {
		a = a % b;
		std::swap(a, b);
	}
	return a;
}
long long lcm(long long a, long long b) { return (a / gcd(a, b)) * b; }
template <typename T>
constexpr T mod(T a, T b) {
	return ((a % b) + b) % b;
}
// Finds x, y s.t. ax + by = d = gcd(a, b).
template <typename T>
constexpr void extended_euclid(T a, T b, T& x, T& y, T& d) {
	if(b == 0)
		x = a > 0 ? 1 : -1, y = 0, d = (a < 0 ? -a : a);
	else
		extended_euclid(b, a % b, y, x, d), y -= a / b * x;
}
template <typename T>
std::pair<T, T> crt(T a1, T m1, T a2, T m2) {
	T s, t, d;
	extended_euclid(m1, m2, s, t, d);
	if(a1 % d != a2 % d) { return {-1, 0}; }
	auto r = mod(s * a2 % m2 * m1 + t * a1 % m1 * m2, m1 * m2) / d;
	auto m = m1 / d * m2;
	return {r, m};
}

// solves ab = 1 (mod n), -1 on failure
constexpr long long mod_inverse(long long a, long long n) {
	long long x = 0, y = 0, d = 0;
	extended_euclid(a, n, x, y, d);
	return (d > 1 ? -1 : mod(x, n));
}
// Along long modular inverses of [1..n] mod P in O(n) time.
std::vector<long long> inverses(long long n, long long P) {
	std::vector<long long> I(n + 1, 1);
	for(long long i = 2; i <= n; ++i) I[i] = mod(-(P / i) * I[P % i], P);
	return I;
}
// (a*b)%m
long long mulmod(long long a, long long b, long long m) {
	return __int128(a) * b % m;
	b %= m;
	long long x = 0, y = a % m;
	while(b > 0) {
		if((b & 1) != 0) x = (x + y) % m;
		y = (2 * y) % m, b /= 2;
	}
	return x % m;
}

template <typename T>
struct Poly;

// Finds b^e in O(lg n) time
template <bool no_inverse = false, typename T = long long>
constexpr T pow(T b, long long e, const T& id = T{1}) {
	assert(e >= 0);
	T p = e < 2 ? id : pow<no_inverse>(b * b, e / 2, id);
	return e & 1 ? p * b : p;
}
template <typename T, typename M>
constexpr T pow(T b, long long e, const T& id, const M& m) {
	// if(e < 0) return 1/pow(b, -e, id, m);
	T p = e < 2 ? id : pow(m(b, b), e / 2, id, m);
	return e & 1 ? m(p, b) : p;
}
// Alias to avoid issues with the std pow.
template <typename T>
constexpr T binpow(T b, long long e, T id = T(1)) {
	return ::pow<T>(b, e, id);
}
// Finds b^e % m in O(lg n) time, ensure that b < m to avoid overflow!
template <typename T = long long>
constexpr T powmod(T b, long long e, T m) {
	T p = e < 2 ? T{1} : powmod((b * b) % m, e / 2, m);
	return e & 1 ? p * b % m : p;
}
// Solve ax + by = c, returns false on failure.
bool linear_diophantine(long long a, long long b, long long c, long long& x, long long& y) {
	long long d = gcd(a, b);
	if((c % d) != 0) return false;
	x = c / d * mod_inverse(a / d, b / d);
	y = (c - a * x) / b;
	return true;
}
