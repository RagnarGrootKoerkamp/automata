#pragma once
#include "./numbertheory.cpp"

#include <ostream>

template <auto P, typename B = long long> // prime, primitive root
struct Field {
	constexpr static B p = P;
	using T              = Field;
	using B_             = B; // NOLINT

  private:
	B x;

  public:
	constexpr Field(B x = 0) : x(x % p) {} // NOLINT

	constexpr T operator+(T r) const { return (x + r.x) % p; }
	constexpr T& operator+=(T r) { return x = (x + r.x) % p, *this; }
	constexpr T operator-() const { return (-x) % p; }
	constexpr T operator-(T r) const { return (x - r.x) % p; }
	constexpr T& operator-=(T r) { return x = (x + p - r.x) % p, *this; }
	constexpr T operator*(T r) const { return (x % p * r.x) % p; }
	constexpr T& operator*=(T r) { return x = (x % p * r.x) % p, *this; }
	constexpr T operator/(T r) const { return (*this) * r.inv(); }
	constexpr T& operator/=(T r) { return (*this) *= r.inv(); }
	constexpr bool operator==(T r) const { return (x - r.x) % p == 0; }

	constexpr B operator%(B r) const {
		assert(P % r == 0);
		return val() % r;
	}

	constexpr friend T operator+(T l, B r) { return (l.x + r) % p; }
	constexpr friend T operator+(B l, T r) { return (l + r.x) % p; }
	constexpr T& operator+=(B r) { return x = (x + r) % p, *this; }
	constexpr friend T operator-(T l, B r) { return (l.x - r + p) % p; }
	constexpr friend T operator-(B l, T r) { return (l - r.x + p) % p; }
	constexpr T& operator-=(B r) { return x = (x + p - r) % p, *this; }
	constexpr friend T operator*(T l, B r) { return (l.x * (r % p)) % p; }
	constexpr friend T operator*(B l, T r) { return ((l % p) * r.x) % p; }
	constexpr T& operator*=(B r) { return x = (x * r) % p, *this; }
	constexpr friend T operator/(T l, B r) { return l * T{r}.inv(); }
	constexpr friend T operator/(B l, T r) { return l * r.inv(); }
	constexpr T& operator/=(B r) { return (*this) *= T{r}.inv(); }
	constexpr bool operator==(B r) const { return (x - r) % p == 0; }

	constexpr T& operator++() {
		++x;
		return *this;
	}

	B val() const { return ((x % p) + p) % p; }
	B val_fast_unsafe() const { return x; }

	constexpr T inv() const {
		// static_assert(IsPrime(p));
		assert(x != 0);
		auto i = mod_inverse(x, p);
		assert(i != -1);
		return i;
	} // pow(x, p-2) is 6x slower
	static T root(B k) {
		assert((p - 1) % k == 0);
		constexpr auto w = primitive_root();
		if(k > 0) return pow(w, (p - 1) / k);
		return pow(w, (p - 1) / -k).inv();
	}
	static constexpr T primitive_root() {
		//static_assert(IsPrime(p));
		T w = 2;
		while(true) {
			auto r  = p - 1;
			bool ok = true;
			for(long long q = 2; q * q <= r; ++q)
				if(r % q == 0) {
					while(r % q == 0) r /= q;
					if(pow(w, (p - 1) / q).x == 1) {
						ok = false;
						break;
					}
				}
			if(r > 1 && pow(w, (p - 1) / r).x == 1) ok = false;
			if(ok) return w.x;
			w += 1;
		}
	}
	[[nodiscard]] bool zero() const { return x == 0LL; }
	// operator B(){return x;}
	constexpr T sqrt() const {
		auto y = sqrtp(x, p);
		assert(y >= 0);
		return {y};
	}
	friend T sqrt(T t) { return t.sqrt(); }
	friend std::ostream& operator<<(std::ostream& o, T t) { return o << (t.x + p) % p; }

	struct comp {
		bool operator()(const T& l, const T& r) const { return l.val() < r.val(); }
	};

	template <B Q>
	constexpr Field<P * Q, B> Multiply() {
		return x * Q;
	}

	template <auto D>
	constexpr Field<P / D, B> Div() {
		static_assert(P % D == 0);
		assert(x % D == 0);
		return x / D;
	}

	template <auto, typename>
	friend struct Field;

	template <auto, auto, typename, typename>
	friend struct FieldFrac;

	// For sorting in std::set.
	friend bool operator<(const Field& l, const Field& r) { return l.val() < r.val(); }
};

template <auto P, typename B = long long> // prime, primitive root
auto abs(Field<P, B> f) {
	return f.val();
}

using F1 = Field<1'107'296'257>; // k * 2^25 + 1
using F2 = Field<2'013'265'921>; // k * 2^27 + 1
using F3 = Field<2'281'701'377>; // k * 2^27 + 1

// Specializations of Pow for Field.
template <auto M, typename B = long long, auto M2 = Phi(M), typename B2 = long long>
constexpr Field<M> pow(Field<M, B> b, Field<M2, B2> e) {
	static_assert(M2 % Phi(M) == 0);
	return pow(b, e.x);
}

// TODO: Make this work for non-coprime A and B.
template <typename T = long long, T P, T Q>
constexpr auto crt(Field<P, T> a, Field<Q, T> b) {
	using A = Field<P, T>;
	using B = Field<Q, T>;
	static_assert(Gcd(A::p, B::p) == 1);
	return (a / B::p).template Multiply<B::p>() + (b / A::p).template Multiply<A::p>();
}
template <typename T = long long, T P, T Q>
constexpr auto operator&(Field<P, T> a, Field<Q, T> b) {
	return crt(a, b);
}
