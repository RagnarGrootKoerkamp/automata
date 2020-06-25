#pragma once
#include "field.cpp"
#include "numbertheory.cpp"
#include "printing.cpp"

#include <cassert>
#include <future>
#include <type_traits>
#include <utility>
#include <vector>

void brinc(int& x, int k) {
	int i = k - 1, s = 1 << i;
	x ^= s;
	if((x & s) != s) {
		--i;
		s >>= 1;
		while(i >= 0 && ((x & s) == s)) x = x & ~s, --i, s >>= 1;
		if(i >= 0) x |= s;
	}
}

template <typename T>
void fft(std::vector<T>& A, int p, bool inv = false) {
	int N = 1 << p;
	for(int i = 0, r = 0; i < N; ++i, brinc(r, p))
		if(i < r) std::swap(A[i], A[r]);
	for(int m = 2; m <= N; m <<= 1) {
		T w, w_m = T::root(inv ? -m : m);
		for(int k = 0; k < N; k += m) {
			w = T{1};
			for(int j = 0; j < m / 2; ++j) {
				T t              = w * A[k + j + m / 2];
				A[k + j + m / 2] = A[k + j] - t;
				A[k + j]         = A[k + j] + t;
				w                = w * w_m;
			}
		}
	}
	if(inv) {
		T inverse = T(N).inv();
		for(auto& x : A) x = x * inverse;
	}
}

// convolution leaves A and B in frequency domain state.
template <typename T>
std::vector<T> convolution(std::vector<T>& A, std::vector<T>& B) {
	int s = A.size() + B.size() - 1;
	if(s == 1) return {A[0] * B[0]};
	int q = 32 - __builtin_clz(s - 1), N = 1 << q;
	A.resize(N, {});
	B.resize(N, {});
	fft(A, q, false);
	fft(B, q, false);
	std::vector<T> C(N);
	for(int i = 0; i < N; ++i) C[i] = A[i] * B[i];
	fft(C, q, true);
	C.resize(s);
	return C;
}

template <typename T>
void pow(std::vector<T>& A, long long e) {
	if(A.size() == 1) return A[0] = pow(A[0], e), void();
	int s = e * (A.size() - 1) + 1, q = 32 - __builtin_clz(s - 1), N = 1 << q;
	A.resize(N, {});
	fft(A, q, false);
	for(auto& x : A) x = pow(x, e);
	fft(A, q, true);
	A.resize(s);
}

// Solves x = ai mod mi. x is unique modulo lcm mi.
// Returns {0, -1} on failure, {x, lcm mi} otherwise.
template <typename T>
std::pair<T, T> crt(const std::vector<T>& a, const std::vector<T>& m) {
	std::pair<T, T> res = {a[0], m[0]};
	for(int i = 1; i < a.size(); ++i) {
		res = crt<T>(res.first, res.second, mod(a[i], m[i]), m[i]);
		if(res.second == -1) break;
	}
	return res;
}

template <typename T>
struct Poly;

template <auto M, typename B>
Poly<Field<M, B>> convolution(const Poly<Field<M, B>>& a, const Poly<Field<M, B>>& b) {
	static_assert(M == 2);
	assert((M - 1) * (M - 1) * std::min(a.size(), b.size()) <= F1::p);

	std::vector<F1> A2(a.size()), B2(b.size());
	// Make sure to use positive representations for all modulo values.
	for(int i = 0; i < a.size(); ++i) A2[i] = a[i].val();
	for(int i = 0; i < b.size(); ++i) B2[i] = b[i].val();
	auto C1 = convolution(A2, B2);

	Poly<Field<M, B>> v(a.size() + b.size() - 1);
	for(int i = 0; i < v.size(); ++i) v[i] = C1[i].val() % M;
	return v;
}

template <long long M, typename B>
Poly<Field<M, B>> fftpow(const Poly<Field<M, B>>& A, long long e) {
	static_assert(M == 2);
	// TODO: For larger exponents, we need more fields / in between steps.
	assert(e <= 2);
	assert((M - 1) * (M - 1) * A.size() <= F1::p);

	std::vector<F1> A2(A.size());
	for(int i = 0; i < A.size(); ++i) A2[i] = A[i].val();
	pow(A2, e);

	Poly<Field<M, B>> v(e * (A.size() - 1) + 1);
	for(int i = 0; i < v.size(); ++i) v[i] = A2[i];

	return v;
}
