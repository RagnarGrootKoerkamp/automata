This program find 2-automata on at most 6 vertices satisfying the following conditions:

- Leading zeros invariant, e.g. outgoing `0` edges point to vertices of the same label.
- Minimal: there exists no smaller automaton with the same power series.
- The power series has `a_0 = 0`, `a_1 = 1` and `a_2 = 1`, i.e. it starts with `\sigma = x + x^2 + O(x^3)`.
- The power series has order `1`, `2`, `4`, or `8` modulo `x^{1025}`.

## Implementation notes

The search over automata works as follows:

- Iterate over the number of vertices `1\leq n \leq 6`. Vertices are numbered `0` to `n-1`. Vertex `0` is the start, and the `1`-edge out of `0` always points to vertex `n-1`.
- Iterate over `1\leq z \leq n-1`, the number of vertices with label `0`.
- Vertices 0 to `z-1` get label `0`, vertices `z` to `n-1` get label `1`.
- Iterate over all possible combinations of `0`- and `1`-edges out of all vertices, where `0` edges must go to other vertices of the same label.
- To break symmetry, the `0`-edge out of vertex `0` must go to either vertex `0` or vertex `1`, and the `0`-edge out of `n-1` must go to either `n-2` or `n-1`.
- For each automaton, compute the power series `\sigma` up to degree `1024`. Then compute `\sigma^8 \mod x^{1025}`. If `\sigma^8 = x`, we check whether `\sigma^1 = x`, `\sigma^2=x`, and `\sigma^4=x` to determine the order.

## Output format

Each automaton in the output first lists the labels of the vertices, followed by a description of the edges. The `i`th line contains the destination of the `0`-edge followed by the destination of the `1`-edges.

For example:
```
Vertex labels: 0 0 1 1 1
Edges:
0 4
1 2
4 2
3 4
3 1
```
This automaton has 5 vertices, where vertices `0` and `1` have label `0` and vertices `2`, `3`, and `4` have label `1`.
The `0`-edge out of `0` goes to `0`, and the `1`-edge out of `0` goes to `4`.

## Output data

* [results\_6.txt](results_6.txt): degree up to 6, order up to 8, modulo `x^{2^{10}+1}=x^{1025}`, and only of the form `x+x^2 + O(x^3)`.
* [results\_6\_order\_8\_degree\_65536.txt](results_6_order_8_degree_65536.txt): degree up to 6, order up to 8, modulo `x^{2^{16}+1}=x^{65537}`, and `\sigma` of the form `x+x^2+O(x^3)` or `x+x^4+O(x^5)`. (Power series are only printed to order `x^{1024}`.)

  It finds the following results. Note that the order 8 power series may be false positives because we only check this modulo `x^{65537}`.
    * Number of automata on 2 vertices of order 2: 1
    * Number of automata on 4 vertices of order 2: 1
    * Number of automata on 4 vertices of order 8: 1
    * Number of automata on 5 vertices of order 2: 2
    * Number of automata on 5 vertices of order 4: 1
    * Number of automata on 5 vertices of order 8: 3
    * Number of automata on 6 vertices of order 2: 7
    * Number of automata on 6 vertices of order 4: 1
    * Number of automata on 6 vertices of order 8: 3

* [results\_6\_order\_8\_degree\_2_18.txt](results_6_order_8_degree_2_18.txt): degree up to 6, order up to 8, modulo `x^{2^{18}+1}`, and `\sigma` of the form `x+x^2+O(x^3)` or `x+x^4+O(x^5)`. (Power series are only printed to order `x^{1024}`.)

  It finds the following results. This shows that there are no order 8 power series of the form `x+x^2` or `x+x^4`.
    * Number of automata on 2 vertices of order 2: 1
    * Number of automata on 4 vertices of order 2: 1
    * Number of automata on 4 vertices of order 8: 1
    * Number of automata on 5 vertices of order 2: 2
    * Number of automata on 5 vertices of order 4: 1
    * Number of automata on 6 vertices of order 2: 7
    * Number of automata on 6 vertices of order 4: 1

* [results_6_order_4_degree_2_18_any_first_coefficient.txt](results_6_order_4_degree_2_18_any_first_coefficient.txt): degree up to 6, order up to 4, modulo `x^{2^{18}+1}`, and `\sigma` of any. (Power series are only printed to order `x^{1024}`.)

  It finds the following results. Some of these will definitely be false positives.
    * Number of automata on 2 vertices of order 2: 1
    * Number of automata on 4 vertices of order 1: 1
    * Number of automata on 4 vertices of order 2: 2
    * Number of automata on 5 vertices of order 2: 4
    * Number of automata on 5 vertices of order 4: 1
    * Number of automata on 6 vertices of order 2: 11
    * Number of automata on 6 vertices of order 4: 1

  These same results hold modulo `x^{2^20+1}` as well, see [results_6_order_4_degree_2_20_any_first_coefficient.txt](results_6_order_4_degree_2_20_any_first_coefficient.txt)




## Some open questions

* To which degree `d` do we need to check that `\sigma^k \equiv x \mod x^d` to be sure equality holds as well.
* Can we find an automaton for `\sigma^2`?
* More general, is there an automaton to automaton mapping that gives the automaton corresponding to `\sigma^2`?
