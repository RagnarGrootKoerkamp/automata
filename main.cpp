#include <array>
#include <bit>
#include <cassert>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <vector>

// An small Field class.
template <auto P, typename B> // prime, underlying integer type.
struct Field {
	constexpr static B p = P;
	using T              = Field;

  private:
	B x;

  public:
	constexpr Field(B x = 0) : x{x % p} {} // NOLINT

	constexpr T operator+(T r) const { return (x + r.x) % p; }
	constexpr T& operator+=(T r) { return x = (x + r.x) % p, *this; }
	constexpr T operator*(T r) const { return (x % p * r.x) % p; }
	constexpr T& operator*=(T r) { return x = (x % p * r.x) % p, *this; }
	constexpr bool operator==(T r) const { return (x - r.x) % p == 0; }
	constexpr bool operator!=(T r) const { return !(*this == r); }

	[[nodiscard]] B val() const { return ((x % p) + p) % p; }

	friend std::ostream& operator<<(std::ostream& o, T t) { return o << t.val(); }

	// Comparison operator so we can use this in std::set.
	bool operator<(const T& r) const { return val() < r.val(); };
};

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
		for(int i = 0; i < size(); ++i)
			for(int j = 0; j < r.size(); ++j) a[i + j] += operator[](i) * r[j];
		return a;
	}
	Poly& operator*=(const Poly& r) { return *this = *this * r; }

	// Note: square means substitution into itself!
	[[nodiscard]] Poly square() const {
		Poly a{0}, pp{1};
		for(const auto& x : *this) {
			if(x != 0) a += x * pp;
			pp *= *this;
			if(pp.size() > size()) pp.resize(size());
		}
		return a;
	}

	// Check whether this polynomial is the identity under substitution: x.
	[[nodiscard]] bool is_identity() const {
		if(size() < 2) return false;
		if(get(0) != 0) return false;
		if(get(1) != 1) return false;
		for(int i = 2; i < size(); ++i)
			if(get(i) != 0) return false;
		return true;
	}

	// The order under substitution of this polynomial, modulo x^1025
	[[nodiscard]] int order() const {
		Poly p, ps, pss, psss;
		for(int d : {4, 8, 16, 32, 64, 128, 256, 512, 1024}) {
			assert(d <= size());
			p = *this;
			p.resize(d);
			ps   = p.square();
			pss  = ps.square();
			psss = pss.square();
			if(not pss.is_identity()) return -1;
		}
		if(p.is_identity()) return 1;
		if(ps.is_identity()) return 2;
		if(pss.is_identity()) return 4;
		if(psss.is_identity()) return 8;
		assert(false);
	}

	// Print the polynomial.
	friend std::ostream& operator<<(std::ostream& o, const Poly& p) {
		bool f = true;
		for(int i = 0; i < p.size(); ++i) {
			if(p[i] == 0) continue;
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
			o << c;
			if(i == 1) o << "*x";
			if(i > 1) o << "*x^" << i;
		}
		return o;
	}

  private:
	[[nodiscard]] const T& get(int i) const { return (*this)[i]; }
	T& get(int i) { return (*this)[i]; }
};

////////////////////////

using F           = Field<2, int>;
using PowerSeries = Poly<F>;

// Automaton on N vertices, where each vertex has 2 outgoing edges: 0 and 1.
// 0 edges must go to the same value.
// start vertex (vertex 0) has value 0.
template <int N>
struct Automaton {
	std::array<F, N> vertex_label;
	std::array<std::array<int, 2>, N> edges{};

	static constexpr int start = 0;

	// Check some conditions:
	// - a_0 = 0
	// - a_1 = 1
	// - a_2 = 1
	// - all vertices are reachable
	// - each vertex with outgoing edges can reach the opposite label.
	[[nodiscard]] bool check() const {
		if(coefficient(0) != 0) return false;
		if(coefficient(1) != 1) return false;
		if(coefficient(2) != 1) return false;

		// Check that vertices are distinct.
		for(int i = 0; i < N; ++i)
			for(int j = 0; j < i; ++j)
				if(vertex_label[i] == vertex_label[j] and edges[i] == edges[j]) return false;

		// Check that all vertices are reachable.
		std::array<bool, N> done{false};
		for(int i = 0; i < (1 << N); ++i) {
			int n   = i;
			int v   = start;
			done[v] = true;
			while(n != 0) {
				v       = step(v, n % 2);
				done[v] = true;
				n /= 2;
			}
		}
		for(int i = 0; i < N; ++i)
			if(not done[i]) return false;

		// Check that each vertex that has outgoing edges can reach vertices of the opposite label.
		for(int i = 0; i < N; ++i) {
			if(edges[i][0] == i and edges[i][1] == i) continue;

			std::array<bool, N> done{};
			std::function<bool(int)> dfs;
			dfs = [&](int u) {
				if(done[u]) return false;
				done[u] = true;
				if(vertex_label[u] != vertex_label[i]) return true;
				for(auto v : edges[u])
					if(dfs(v)) return true;
				return false;
			};

			if(not dfs(i)) return false;
		}

		return true;
	}

	// Walking the automaton.
	[[nodiscard]] int step(int v, F bit) const { return edges[v][bit.val()]; }

	// Get a_n, the coefficient of x^n.
	[[nodiscard]] F coefficient(int n) const {
		int v = start;
		while(n != 0) {
			v = step(v, n % 2);
			n /= 2;
		}
		return vertex_label[v];
	}

	// Return the power series for this automaton, up to and including x^degree.
	[[nodiscard]] PowerSeries power_series(int degree) const {
		// Make sure we take enough steps to cover all vertices.
		PowerSeries p(degree + 1);
		for(int i = 0; i <= degree; ++i) p[i] = coefficient(i);
		return p;
	}

	friend std::ostream& operator<<(std::ostream& o, const Automaton& a) {
		o << "Vertex labels: ";
		for(auto x : a.vertex_label) o << x << " ";
		o << std::endl;
		o << "Edges:\n";
		for(auto& es : a.edges) {
			for(auto e : es) o << e << " ";
			o << std::endl;
		}
		return o;
	}
};

std::set<PowerSeries> seen;

// Count the number of suitable automata.
template <int N>
void count() {
	int num_automata    = 0;
	int num_ok_automata = 0;
	std::map<int, int> count_per_order;

	Automaton<N> automaton;

	// Loop over the number of 1 labels.
	// Note: 0 has label 0, and at least one vertex has label 1, because coefficient(1) = 1.
	for(int zeros = 1; zeros < N; ++zeros) {
		// 0..(zeros-1) have label 0,
		// zeros..(N-1) have label 1.
		for(int i = 0; i < N; ++i) automaton.vertex_label[i] = i < zeros ? 0 : 1;

		// Now iterate over all possible assignments of the edges:
		// The 0 edge goes to another vertex of the same label.
		// The 1 edge goes to an arbitrary vertex.

		// u is the vertex we're currently iterating over.
		std::function<void(int)> dfs;
		dfs = [&](int u) {
			if(u == N) {
				// Check whether the current automaton works.
				++num_automata;

				if(not automaton.check()) return;
				++num_ok_automata;

				auto p = automaton.power_series(1024);
				if(seen.contains(p)) return;
				auto order = p.order();
				if(order <= 0) return;
				auto [it, inserted] = seen.insert(p);
				if(not inserted) return;
				++count_per_order[order];
				std::cout << "Found solution of order " << order << std::endl;
				std::cout << automaton << std::endl;
				std::cout << "PowerSeries: " << p << std::endl << std::endl << std::endl;
				return;
			}

			// Loop over 0-edge to vertex of same label.
			auto& v0 = automaton.edges[u][0];
			auto& v1 = automaton.edges[u][1];

			// 0-edge goes to any other vertex with the same label.
			// For the start vertex, always use 0 or 1 to break symmetry.
			// For the final vertex, always use N-2 or N-1 to break symmetry.
			for(v0 = (u == N - 1 ? N - 2 : 0); v0 < (u == 0 ? 2 : N); ++v0) {
				if(automaton.vertex_label[v0] != automaton.vertex_label[u]) continue;

				// Loop over 1-edge to arbitrary vertex.
				// For the start vertex, we force the edge to go to vertex N-1 to break symmetry.
				if(u == 0)
					for(v1 = N - 1; v1 < N; ++v1) dfs(u + 1);
				else
					for(v1 = 0; v1 < N; ++v1) dfs(u + 1);
			}
		};

		dfs(0);
	}

	std::cout << "Number of automata on " << N
	          << " vertices:                        " << num_automata << std::endl;
	std::cout << "Number of automata on " << N
	          << " vertices post check:             " << num_ok_automata << std::endl;
	for(auto [order, cnt] : count_per_order)
		std::cout << "Number of automata on " << N << " vertices of order " << order << ": " << cnt
		          << std::endl;
	std::cout << std::endl;

	std::cout << "========================================================" << std::endl << std::endl;
}

int main() {
	count<1>();
	count<2>();
	count<3>();
	count<4>();
	count<5>();
	count<6>();
}
