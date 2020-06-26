#include "lib/fft.cpp"
#include "lib/field.cpp"
#include "lib/poly.cpp"

#include <array>
#include <bit>
#include <cassert>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <vector>

////////////////////////

using F           = Field<2, char>;
using PowerSeries = Poly<F>;

constexpr int FAST_DEGREE = 1 << 8;
constexpr int FULL_DEGREE = 1 << 20;
// max order = 2^MAX_LOG_ORDER
constexpr int MAX_LOG_ORDER = 2;

// Operations on compositional power series.

// Using n n*lg n FFT multiplications.
[[nodiscard]] PowerSeries square_fft(const PowerSeries& p) {
	assert(p[0] == 0);
	assert(p[1] == 1);

	PowerSeries pp{1};
	PowerSeries ans(p.size());
	for(const auto& ai : p) {
		if(ai != 0) ans += ai * pp;
		pp *= p;
		if(pp.size() > p.size()) pp.resize(p.size());
	}
	return ans;
}

// Multiply by largest 2^a <= i.
[[nodiscard]] PowerSeries square_nlogn(const PowerSeries& p) {
	// return square_fft(p);
	assert(p[0] == 0);
	assert(p[1] == 1);

	// p\circ p = p + a_2 p^2 + a_3 p^3 + ...
	// When k = 2^a+b, with b < 2^a, we have p^k = p^(2^a) * p^b, and p^(2^a) = sum_i a_i x^(i*2^a)
	std::vector<PowerSeries> powers(p.size());
	powers[0] = PowerSeries{1};

	PowerSeries ans(p.size());
	for(unsigned i = 1; i < p.size(); ++i) {
		int a = std::bit_floor(i);
		int b = i - a;
		powers[i].resize(p.size());
		for(int j = 0; a * j < p.size(); ++j) {
			if(p[j] == 0) continue;
			assert(p[j] == 1);
			for(int k = 0; k < std::min(powers[b].size(), p.size() - a * j); ++k)
				powers[i][a * j + k] += powers[b][k];
		}
		if(p[i] != 0) ans += powers[i];
	}
	return ans;
}

// Multiply by 2^b|i, or largest 2^a<=i if b=0.
[[nodiscard]] PowerSeries square_fast_loop(const PowerSeries& p) {
	// return square_fft(p);
	assert(p[0] == 0);
	assert(p[1] == 1);

	// p\circ p = p + a_2 p^2 + a_3 p^3 + ...
	// When k = 2^a+b, with b < 2^a, we have p^k = p^(2^a) * p^b, and p^(2^a) = sum_i a_i x^(i*2^a)
	std::vector<PowerSeries> powers(p.size());
	powers[0] = PowerSeries{1};

	PowerSeries ans(p.size());
	for(unsigned i = 1; i < p.size(); i += 2) {
		int z = std::countr_zero(i);
		powers[i].resize(p.size());
		int a = std::bit_floor(i);
		int b = i - a;
		for(int j = 0; a * j < p.size(); ++j) {
			if(p[j] == 0) continue;
			assert(p[j] == 1);
			for(int k = 0; k < std::min(powers[b].size(), p.size() - a * j); ++k)
				powers[i][a * j + k] += powers[b][k];
		}
		if(p[i] != 0) ans += powers[i];
		for(z = 2; i * z < p.size(); z *= 2) {
			if(p[z * i] == 1) {
				for(int j = 0; z * j < p.size(); ++j) ans[j * z] += powers[i][j];
			}
		}
	}
	return ans;
}

// Multiply by 2^b|i, or largest 2^a<=i if b=0.
[[nodiscard]] PowerSeries square_dfs(const PowerSeries& p) {
	// return square_fft(p);
	assert(p[0] == 0);
	assert(p[1] == 1);

	// p\circ p = p + a_2 p^2 + a_3 p^3 + ...
	// When k = 2^a+b, with b < 2^a, we have p^k = p^(2^a) * p^b, and p^(2^a) = sum_i a_i x^(i*2^a)
	PowerSeries ans(p.size());

	std::vector<bool> needed(p.size(), false);
	needed[0] = true;
	for(int i = 1; i < p.size(); ++i) {
		if(p[i] == 0) continue;
		unsigned j = i;
		while(j % 2 == 0 and not needed[j]) {
			needed[j] = true;
			j /= 2;
		}
		while(not needed[j]) {
			assert(j > 0);
			needed[j] = true;
			j         = j - std::bit_floor(j);
		}
	}

	std::function<void(int, const PowerSeries&, int)> dfs;
	std::vector<PowerSeries> cache(std::bit_width(p.size()) + 2);
	// Compute power series for i using given power series for i - bit_floor(i).
	dfs = [&](unsigned i, const PowerSeries& parent, int d = 0) {
		assert(i < p.size());
		if(not needed[i]) return;
		PowerSeries& current = cache[d];
		current.assign(p.size(), 0);
		const int a = std::bit_floor(i);
		const int b = i - a; // corresponds to parent.
		for(int j = 0; a * j < p.size(); ++j) {
			if(p[j] == 0) continue;
			assert(p[j] == 1);
			for(int k = 0; k < std::min(parent.size(), p.size() - a * j); ++k)
				current[a * j + k] += parent[k];
		}
		if(p[i] != 0) ans += current;
		for(int z = 2; i * z < p.size(); z *= 2) {
			if(p[z * i] == 1) {
				for(int j = 0; z * j < p.size(); ++j) ans[j * z] += current[j];
			}
		}

		for(int aa = 2 * a; i + aa < p.size(); aa *= 2) dfs(i + aa, current, d + 1);
	};

	dfs(1, PowerSeries{1}, 0);

	return ans;
}
[[nodiscard]] PowerSeries square(const PowerSeries& p) { return square_dfs(p); }

// Check whether this polynomial is the identity under substitution: x.
[[nodiscard]] bool is_identity(const PowerSeries& p) {
	if(p.size() < 2) return false;
	if(p[0] != 0) return false;
	if(p[1] != 1) return false;
	for(int i = 2; i < p.size(); ++i)
		if(p[i] != 0) return false;
	return true;
}

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

		// Make sure that x^2 or x^4 is the next non-0 coefficient.
		/*
		bool ok = false;
		if(coefficient(2) == 1)
		    ok = true;
		else {
		    ok = coefficient(3) == 0 and coefficient(4) == 1;
		}
		if(not ok) return false;
		*/

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

// The order under substitution of this polynomial, modulo x^1025
template <int N>
[[nodiscard]] int order(const Automaton<N>& a, PowerSeries& p, int max_deg) {
	std::array<PowerSeries, MAX_LOG_ORDER + 1> powers;
	p.resize(2);
	for(int d = 4; d <= max_deg; d *= 2) {
		// if(d > 8 * FAST_DEGREE) std::cerr << "Degree: " << d << std::endl;
		while(p.size() < d) p.push_back(a.coefficient(p.size()));
		powers[0] = p;
		for(int i = 0; i < MAX_LOG_ORDER; ++i) powers[i + 1] = square(powers[i]);
		if(not is_identity(powers.back())) return -1;
	}
	for(int i = 0; i <= MAX_LOG_ORDER; ++i)
		if(is_identity(powers[i])) return 1 << i;
	assert(false);
}

std::set<PowerSeries> seen;

// Count the number of suitable automata.
template <int N>
void count() {
	int num_automata    = 0;
	int num_ok_automata = 0;
	std::map<int, int> count_per_order;

	PowerSeries p;
	p.reserve(FULL_DEGREE);
	p.push_back(0);
	p.push_back(1);
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

				if(order(automaton, p, FAST_DEGREE) <= 0) return;
				if(seen.contains(p)) return;
				auto order = ::order(automaton, p, FULL_DEGREE);
				if(order <= 0) return;
				auto p2 = p;
				p2.resize(FAST_DEGREE);
				auto [it, inserted] = seen.insert(p2);
				// The FAST_DEGREE terms should be unique already.
				assert(inserted);
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
			for(v0 = (u == N - 1 ? N - 2 : 0); v0 < (u == 0 ? 2 : N); ++v0) { // NOLINT
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

	std::cout << "========================================================" << std::endl
	          << std::endl;
}

int main() {
	// count<1>();
	// count<2>();
	// count<3>();
	// count<4>();
	count<5>();
	count<6>();
}
