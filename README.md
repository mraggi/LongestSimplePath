# LongestSimplePath

Finds long simple paths (no repeated vertices) in large weighted digraphs: given arc weights, possibly negative, it searches for a simple path of maximum total weight. With unit weights this is the longest path problem. It also computes an upper bound, so it can sometimes prove that the path it found is optimal. Paths can optionally be required to start at a given vertex.

A previous version of this code won the Oracle MDC 2015 coding challenge.

## Build

Requirements: a C++20 compiler (tested with GCC 16), CMake ≥ 3.20, and optionally GoogleTest for the tests.

```sh
cmake -S . -B build -G Ninja        # Release (-O3) by default
ninja -C build
./build/lsp_tests                   # unit tests (needs GoogleTest)
./build/lsp_example                 # minimal library example
```

## Usage

```sh
./build/lsp solve graph.el --path-out best.path       # search for 10 s (the default), write the best path
./build/lsp solve graph.el --time 60 --source 0        # search for 60 s; paths must start at vertex 0
./build/lsp ap graph.el                                # only the upper bound and the fast initial path
./build/lsp info graph.el                              # size, strong components, condensation depth
./build/lsp_check graph.el best.path                   # verify a path file
```

The search stops at the time limit (`--time`, default 10 s; `--time 0` means no limit), or earlier if the path reaches the upper bound, which proves it optimal. `--seed S` (default 1) fixes the random choices.

**Input.** By default an edge list with 0-based vertices: a header line `n m`, then one line `u v [w]` per arc (weight 1 if omitted). Other formats, with `--format`:
* `adjlist`: the 1-based adjacency lists of Farnebäck's [Go Game Graphs](https://github.com/GunnarFarneback/GoGameGraphs);
* `snap`: SNAP edge lists with arbitrary vertex ids and `#` comments;
* `dimacs`: `p sp` / `a u v w`.

Parallel arcs are merged (keeping the maximum weight) and loops are ignored.

**Output.** `lsp solve` prints one JSON object. Its fields include:
* `value`, and `check` (`ok` if the returned path was verified);
* `ap_bound`: the upper bound (valid when `ap_certified` and `ap_bound_global` are true);
* `proven_optimal`;
* `trajectory`: the time and value of every improvement;
* timings and counters;
* with `--print-path`, the path itself.

**Methods.** The default is the best method of this program, `--preset newap`: the assignment relaxation followed by the local search (see [How it works](#how-it-works)). Alternatives:

| option | what it is |
|---|---|
| `--preset new` | the local search alone, without the assignment relaxation |
| `--preset 2018` | the 2016–2018 method |
| `--method rdfs`, `greedy`, `gamm` | baselines: randomized DFS, Warnsdorff greedy, the mutation GA of Portugal, Antunes & Rocha (2010) |

Every component can be switched with `--set key=value`; see [`docs/OPTIONS.md`](docs/OPTIONS.md).

### As a library

```cpp
#include "lsp/solver.hpp"

lsp::Digraph g(6, {{0, 1, 3}, {1, 2, 2}, {2, 0, 4}, {2, 3, 1}, {3, 4, 5}, {1, 4, 1}, {4, 5, 2}, {5, 3, -1}});
lsp::SolverOptions opt;
// default-constructed options are the recommended method (preset "newap"), 10 s
opt.time_limit = 1;             // seconds
// opt.source = 0;              // paths must start at vertex 0
lsp::SolveResult res = lsp::solve(g, opt);
// res.best.vertices, res.best.value, res.proven_optimal, res.ap.bound (if res.ap.certified && res.ap.bound_global)
```

Link against the `lsp` CMake target. The full program is [`examples/minimal.cpp`](examples/minimal.cpp).

### Python scripts (`scripts/`)

* `pathcheck.py GRAPH PATHFILE [SOURCE]`: an independent checker sharing no code with the solver (path file: whitespace-separated 0-based vertex ids).
* `verify_ap_bound.py GRAPH [--source s]`: recomputes the upper bound with SciPy, as an independent check of `ap_bound`.
* `hybrid.py GRAPH --time T --seed S --lsp build/lsp [--source s]`: runs `lsp` and then hands its path to the CP-SAT solver as a hint, on graphs with at most 10⁴ vertices (needs `ortools`). On small random digraphs this finds better paths and can prove optimality.
* `cpsat_heuristic.py`: the CP-SAT model used by `hybrid.py`, also usable on its own.

## How it works

The search combines two ingredients that are strong on different kinds of graphs.

1. **Assignment relaxation.** Every vertex chooses a successor, itself (unused) or "end of path", with each vertex chosen at most once. The best such assignment is an upper bound on the weight of every simple path, for arbitrary weights. It is a collection of one path and some cycles. The code solves it with an ε-scaling auction and checks its optimality exactly with an integer dual certificate. It then patches the cycles into the path. On sparse directed graphs without much symmetry, such as random digraphs, this gives a path close to optimal in a fraction of a second. On symmetric graphs (grids, geometric graphs, undirected graphs) almost every vertex sits on a 2-cycle and the relaxation is of little use.
2. **Local search on pseudo-topological orders.** A pseudo-topological order is a linear order of the vertices in which every backward arc stays inside one strongly connected component. Over such an order, a dynamic program finds the best path that uses only forward arcs, in linear time. The search keeps the current path as a subsequence of the order. It moves the vertices of one component that are not on the path into a gap of the path, ordered by the reverse post-order of a depth-first search, and lets the dynamic program find the best detour through them. This is strong exactly where the relaxation is weak: geometric graphs, grids and the Go game graphs.

An iterated local search repeats this until the time limit: it keeps a random piece of the best path, extends it by depth-first search and improves it again.

## References

* M. Raggi, *Finding long simple paths in a weighted digraph using pseudo-topological orderings*, [arXiv:1609.07450](https://arxiv.org/abs/1609.07450) (2016). A new paper describing the current method and an experimental study is in preparation.
* G. Farnebäck, [Go Game Graphs](https://github.com/GunnarFarneback/GoGameGraphs): large directed longest-path instances with published bounds.

Changes and bug fixes are listed in [`CHANGELOG.md`](CHANGELOG.md).

## License

Apache License 2.0 (see `LICENSE`). Research supported in part by PAPIIT IA106316, UNAM.
