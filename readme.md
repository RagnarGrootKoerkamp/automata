This program find 2-automata on at most 6 vertices satisfying the following conditions:

- Leading zeros invariant, e.g. outgoing `0` edges point to vertices of the same label.
- Minimal: there exists no smaller automaton with the same power series.
- The power series has `a_0 = 0`, `a_1 = 1` and `a_2 = 2`, i.e. it starts with `\sigma = x + x^2 + O(x^3)`.
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

The file [results_6.txt](results_6.txt) contains results for up to `6` vertices. It prints all automata with degree at most `8` and their corresponding power series. Each automaton first lists the labels of the vertices, followed by a description of the edges. The `i`th line contains the destination of the `0`-edge followed by the destination of the `1`-edges.

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
