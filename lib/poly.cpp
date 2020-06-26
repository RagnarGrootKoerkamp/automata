#pragma once
#include "fft.cpp"

#include <vector>

// Polynomials.
template <typename T>
struct Poly : public std::vector<T> {
	using std::vector<T>::vector, std::vector<T>::size, std::vector<T>::operator[];
	explicit Poly(std::vector<T>&& v) : std::vector<T>(std::move(v)) {}
	explicit Poly(const std::vector<T>& v) : std::vector<T>(v) {}
	explicit Poly(const T& t) : std::vector<T>{t} {}

	static const Poly x;

	Poly<T>& resize(int n) {
		std::vector<T>::resize(n);
		return *this;
	}

	// Remove trailing zeros.
	Poly<T>& simplify() {
		while(!this->empty() && this->back() == 0) this->pop_back();
		return *this;
	}

	// Addition.
	Poly<T> operator+=(const Poly<T>& r) {
		auto& l = *this;
		if(r.size() > size()) l.resize(r.size());
		for(int i = 0; i < r.size(); ++i) l[i] += r[i];
		return l;
	}
	friend Poly<T> operator+(Poly<T> l, const Poly<T>& r) {
		if(r.size() > l.size()) l.resize(r.size());
		for(int i = 0; i < r.size(); ++i) l[i] += r[i];
		return l;
	}

	// Linear time multiplication by a constant.
	friend Poly<T> operator*(const T& a, Poly<T> r) {
		for(auto& x : r) x *= a;
		return r;
	}
	friend Poly<T> operator*(Poly<T> r, const T& a) { return a * r; }

	// Quadratic time multiplication.
	// TODO(ragnar): Use Karatsuba or FFT so we can check higher degree.
	Poly operator*(const Poly& r) const {
		if(this->empty() or r.empty()) return {};
		Poly a(this->size() + r.size() - 1);
		// For large polynomials, use FFT.
		constexpr int S = 200;
		if(this->size() > S and r.size() > S) return convolution(*this, r);

		// For small polynomials, use quadratic multiplication.
		for(int i = 0; i < size(); ++i)
			for(int j = 0; j < r.size(); ++j) a[i + j] += operator[](i) * r[j];
		return a;
	}
	Poly& operator*=(const Poly& r) { return *this = *this * r; }

	// Print the polynomial.
	friend std::ostream& operator<<(std::ostream& o, const Poly& p) {
		bool f  = true;
		int cnt = 0;
		for(int i = 0; i < p.size(); ++i) {
			if(p[i] == 0) continue;
			if(++cnt > 500) break;
			auto c = p[i];
			if(f)
				f = false;
			else {
				o << " ";
				if constexpr(std::is_integral_v<T>) {
					if(c > 0)
						o << "+";
					else {
						o << "-";
						c = -c;
					}
				} else {
					o << "+";
				}
				o << " ";
			}
			if(c != 1) o << c << "*";
			if(i == 1) o << "x";
			if(i > 1) o << "x^" << i;
		}
		return o;
	}

  private:
	[[nodiscard]] const T& get(int i) const { return (*this)[i]; }
	T& get(int i) { return (*this)[i]; }
};

