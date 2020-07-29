#include "lib/field.cpp"
#include "lib/poly.cpp"

#include <array>
#include <bit>
#include <bitset>
#include <cassert>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

////////////////////////

// Must be a power of 2. The max order to consider.
constexpr unsigned MAX_ORDER = 4;
constexpr int MAX_LOG_ORDER  = std::countr_zero(MAX_ORDER);

constexpr int FULL_DEGREE = (1 << 14);

////////////////////////
static_assert(std::has_single_bit(MAX_ORDER));

using F           = Field<2, char>;
using PowerSeries = Poly<F>;

// Compute the compositional square of p.
// In general p^(2^k) is easy to compute over F_2: we get b_(2^k*i) = a_i and all other b_i are 0.
// We compute p^i for all odd i by looking at the maximal a with 2^a<=i and compute
// p^i = (p^(2^a) * p^(i-2^a)).
// All even p^(i * 2^b) are computed as (p^i)^(2^b).
// Multiply by 2^b|i, or largest 2^a<=i if b=0.
[[nodiscard]] PowerSeries square(const PowerSeries& p) {
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
	// Cache is only to reuse allocated memory.
	std::vector<PowerSeries> cache(std::bit_width(p.size()) + 2);
	// Compute power series for i using given power series for i - bit_floor(i).
	dfs = [&](unsigned i, const PowerSeries& parent, int d = 0) {
		assert(i < p.size());
		if(not needed[i]) return;
		PowerSeries& current = cache[d];
		current.assign(p.size(), 0);
		const int a = std::bit_floor(i);
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

	// See https://en.wikipedia.org/wiki/DFA_minimization#Hopcroft's_algorithm.
	[[nodiscard]] bool is_minimal() const {
		constexpr auto NN  = 1 << N;
		unsigned accepting = 0, rejecting = 0;
		for(int i = 0; i < N; ++i)
			if(vertex_label[i] == 1)
				accepting |= 1 << i;
			else
				rejecting |= 1 << i;
		std::vector<unsigned> P{accepting, rejecting}, W{accepting, rejecting};
		std::bitset<NN> doneW{}, doneP{};
		P.reserve(NN);
		doneP[accepting] = doneP[rejecting] = true;
		while(not W.empty()) {
			auto A = W.back();
			W.pop_back();
			if(doneW[A]) continue;
			doneW[A] = true;
			for(auto c : {0, 1}) {
				unsigned X = 0;
				for(int x = 0; x < N; ++x)
					if(A & (1 << edges[x][c])) X |= 1 << x;
				for(auto& Y : P) {
					const auto I = Y & X;
					const auto J = Y & (~X);
					if(I and J) {
						// Replace/add to W.
						if(auto it = find(begin(W), end(W), Y); it != W.end()) {
							*it = I;
							W.push_back(J);
						} else {
							if(std::popcount(I) <= std::popcount(J))
								W.push_back(I);
							else
								W.push_back(J);
						}
						// Replace in P.
						doneP[Y] = false;
						doneP[I] = true;
						doneP[J] = true;
						Y        = I;
						P.push_back(J);
					}
				}
			}
		}
		return P.size() == N;
	}

	friend bool equivalent(const Automaton& l, const Automaton& r) {
		if(l.vertex_label != r.vertex_label) return false;
		int Z = 0;
		for(auto x : l.vertex_label)
			if(x == 0) ++Z;
		assert(Z > 0 and Z < N);
		// loop over permutations of [1..Z) and [Z..N-1)
		std::vector<int> p(N);
		std::iota(begin(p), end(p), 0);
		do {
			do {
				// Check if l == r under the current permutation.
				for(int u = 0; u < N; ++u)
					for(auto c : {0, 1})
						if(p[l.edges[u][c]] != r.edges[p[u]][c]) goto bad;
				return true;
			bad:;
			} while(std::next_permutation(begin(p) + Z, end(p) - 1));
		} while(std::next_permutation(begin(p) + 1, begin(p) + Z));
		return false;
	}
};

// The order under substitution of this polynomial, modulo x^1025
template <int N>
[[nodiscard]] int order(const Automaton<N>& a, PowerSeries& p, int max_deg) {
	std::array<PowerSeries, MAX_LOG_ORDER + 1> powers;
	p.resize(2);
	for(int d = 4; d <= max_deg; d *= 2) {
		while(p.size() <= d) p.push_back(a.coefficient(p.size()));
		powers[0] = p;
		for(int i = 0; i < MAX_LOG_ORDER; ++i) powers[i + 1] = square(powers[i]);
		if(not is_identity(powers.back())) return -1;
	}
	for(int i = 0; i <= MAX_LOG_ORDER; ++i)
		if(is_identity(powers[i])) return 1 << i;
	assert(false);
}

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

	std::map<PowerSeries, Automaton<N>> seen;

	auto print_powerseries = [](PowerSeries p) {
		std::cout << "sigma^" << 1 << ": " << p << std::endl;
		for(int i = 1; i <= 4; ++i) {
			p = square(p);
			std::cout << "sigma^" << (1 << i) << ": " << p << std::endl;
		}
	};

	// MAIN WORK IS HERE.
	// Check whether the current automaton works.
	auto test_automaton = [&] {
		++num_automata;

		// Some basic checks (sigma = 0 + t + O(t^2); all states are reachable.)
		if(not automaton.check()) return;
		// Make sure the automaton is minimal.
		if(not automaton.is_minimal()) return;
		++num_ok_automata;

		// Compute the order of the automaton, modulo x^FULL_DEGREE.
		auto order = ::order(automaton, p, FULL_DEGREE);
		// If order modulo x^FULL_DEGREE isn't finite, return.
		if(order <= 0) return;

		auto [it, inserted] = seen.emplace(p, automaton);
		if(not inserted) {
			// If the power series already exists, make sure that the corresponding
			// automaton is equivalent (i.e. the solution is unique up to isomorphisms).
			assert(equivalent(automaton, it->second));
			return;
		}

		// We found a new solution!
		++count_per_order[order];
		std::cout << "Candidate automaton of order " << order << std::endl;
		std::cout << automaton << std::endl;
		print_powerseries(p);
		std::cout << std::endl << std::endl;
	};

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
				test_automaton();
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
				// For the start vertex, we force the edge to go to vertex N-1 to break
				// symmetry.
				if(u == 0)
					for(v1 = N - 1; v1 < N; ++v1) dfs(u + 1);
				else
					for(v1 = 0; v1 < N; ++v1) dfs(u + 1);
			}
		};

		dfs(0);
	}

	// std::cout << "Number of automata on " << N << " vertices:                        " <<
	// num_automata << std::endl; std::cout << "Number of automata on " << N << " vertices post
	// check:             " << num_ok_automata << std::endl;
	for(auto [order, cnt] : count_per_order)
		std::cout << "Number of automata on " << N << " vertices of order " << order << ": " << cnt
		          << std::endl;
	std::cout << std::endl;
	std::cout << "========================================================" << std::endl
	          << std::endl;
}

int main() {
	count<2>();
	count<3>();
	count<4>();
	count<5>();
}
