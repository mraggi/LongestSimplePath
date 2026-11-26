# LongestSimplePath

Heuristics and upper bounds for the **maximum-weight simple path** problem in large weighted digraphs: find a path without repeated vertices whose arcs have maximum total weight (with unit weights, the longest path). Arc weights may be negative. A rooted variant (paths must start at a given vertex) is supported throughout.

This is the 2026 rewrite (C++20) of the code that won the Oracle MDC 2015 coding challenge and accompanied the preprint [arXiv:1609.07450](https://arxiv.org/abs/1609.07450). The original 2016–2018 code is kept unchanged in [`legacy_2018/`](legacy_2018/).

## What is new compared with the 2018 code

* **Assignment relaxation (upper bound + starting path).** Every vertex chooses a successor, itself (unused) or "end of path". The optimum of this assignment problem is an upper bound on every simple path, for arbitrary weights. It is solved by an ε-scaling auction whose optimality is checked exactly (integer dual certificate). Its cycles are then patched into one path (Karp-style). On directed graphs of expander type, such as random digraphs or graphs with a planted Hamiltonian path, the patched path is close to optimal after a fraction of a second. If a path reaches the bound, the solver stops with a **proof of optimality**.
* **Better pseudo-topological local search.** The free vertices of a strongly connected component are ordered by the reverse post-order of a depth-first search started next to the gap, instead of the 2018 randomized "heuristic sort". This makes every DFS tree path available to the dynamic program, and is by far the most important improvement on geometric graphs, grids and Go game graphs.
* **Anytime iterated local search.** It uses deterministic, work-based stopping rules and runs until a time limit. Heavier arcs are tried first on weighted graphs, every non-trivial component is opened, and the whole graph is searched, not only its largest weak component.
* **Rooted variant** (`--source s`), **several input formats**, an **independent path checker**, unit tests, and baselines (randomized DFS, Warnsdorff greedy, the genetic algorithm of Portugal et al. 2010).

### Bugs of the 2018 code that are fixed here
1. The first vertex of every path produced by the improvement phase was dropped from the returned path, although its arc was counted in the reported weight. The reported value was therefore too high by one arc.
2. A stale parent pointer in the dynamic program could produce **invalid paths** (consecutive vertices not joined by an arc). This happened in all runs on a graph with negative weights, and on some positively weighted random graphs.
3. The incrementally maintained maximum ignored the unchanged prefix of the order in one case.

The new code verifies every returned path (`check` field of the output).

## Build

Requirements: a C++20 compiler (tested with GCC 16), CMake ≥ 3.20, and optionally GoogleTest for the tests.

```sh
cmake -S . -B build -G Ninja        # Release (-O3) by default
ninja -C build
./build/lsp_tests                   # unit tests
./build/lsp_example                 # minimal library example (examples/minimal.cpp)
```

## Quick start

```sh
./build/lsp solve graph.el --preset newap --time 60 --seed 1 --path-out best.path
./build/lsp solve go.el --preset newap --time 60 --source 0      # paths must start at vertex 0
./build/lsp ap graph.el                                          # only the assignment bound and the patched path
./build/lsp info graph.el                                        # size, SCCs, condensation depth
```

The default input format is an edge list with 0-based vertices: a header line `n m`, then one line `u v [w]` per arc (weight 1 if omitted). Also accepted, with `--format`:
* `adjlist`: Farnebäck's GoGameGraphs adjacency lists, 1-based;
* `snap`: SNAP edge lists with arbitrary ids and `#` comments;
* `dimacs`: `p sp` / `a u v w`.

Parallel arcs are merged (keeping the maximum weight) and loops are ignored.

`lsp solve` prints one JSON object. Its fields include:
* `value` and `check` (`ok` if the returned path was verified);
* `trajectory`: the time and value of every improvement;
* `ap_bound`, `ap_certified`, `ap_bound_global`: the assignment bound;
* `proven_optimal`;
* timings and counters;
* with `--print-path`, the path itself.

### Presets and methods

| preset / method | what it is |
|---|---|
| `--preset newap` (**recommended**) | assignment relaxation + patching, then the iterated pseudo-topological local search; stops early when the path reaches the certified bound |
| `--preset new` | the iterated local search alone (no relaxation) |
| `--preset 2018` | a faithful re-implementation of the 2016–2018 method (without its bugs) |
| `--method rdfs`, `greedy`, `gamm` | baselines: randomized DFS, Warnsdorff greedy, the mutation GA of Portugal, Antunes & Rocha (2010) |

Every component can be switched with `--set key=value`; see [`docs/OPTIONS.md`](docs/OPTIONS.md).

### When to use what
The two ingredients are complementary.
* **Relaxation strong:** on sparse directed graphs without much symmetry, few cycles appear and patching them gives near-optimal paths.
* **Relaxation useless:** on symmetric digraphs (every edge in both directions: grids, geometric graphs, undirected graphs) almost every vertex sits on a 2-cycle.
* **Local search strong:** on symmetric graphs and on cycle-rich graphs such as the Go game graphs.

`newap` combines both and was the best of the methods we tested on most of our benchmark instances. A constraint solver such as CP-SAT (circuit model) is still better on small random digraphs (a few thousand vertices). [`scripts/hybrid.py`](scripts/hybrid.py) runs `newap` and then hands its path to CP-SAT as a hint on graphs with at most 10⁴ vertices (needs `ortools`).

## Library

```cpp
#include "lsp/solver.hpp"

lsp::Digraph g(6, {{0, 1, 3}, {1, 2, 2}, {2, 0, 4}, {2, 3, 1}, {3, 4, 5}, {1, 4, 1}, {4, 5, 2}, {5, 3, -1}});
lsp::SolverOptions opt;
opt.apply_preset("newap");
opt.time_limit = 1;             // seconds
// opt.source = 0;              // rooted variant
lsp::SolveResult res = lsp::solve(g, opt);
// res.best.vertices, res.best.value, res.proven_optimal, res.ap.bound (if res.ap.certified && res.ap.bound_global)
```

Link against the `lsp` CMake target. The full program is [`examples/minimal.cpp`](examples/minimal.cpp).

## Scripts (`scripts/`, Python 3)
* `pathcheck.py GRAPH PATHFILE [SOURCE]`: independent checker sharing no code with the solver (path file: whitespace-separated 0-based vertex ids).
* `verify_ap_bound.py GRAPH [--source s]`: recomputes the assignment bound with SciPy (Jonker–Volgenant), as an independent check of `ap_bound`.
* `hybrid.py GRAPH --time T --seed S --lsp build/lsp [--source s]`: `newap` followed by CP-SAT with a solution hint (`cpsat_polish.py`).
* `cpsat_heuristic.py`: the CP-SAT circuit model used as an anytime heuristic or exact solver.

## Layout
* `include/lsp`, `src`: the library.
* `apps`: the `lsp` and `lsp_check` command-line tools.
* `tests`: GoogleTest unit tests (including brute-force checks of the dynamic program and of the assignment relaxation on small graphs).
* `examples`, `scripts`, `docs`.
* `legacy_2018/`: the original code, unchanged.

## References
* M. Raggi, *Finding long simple paths in a weighted digraph using pseudo-topological orderings*, [arXiv:1609.07450](https://arxiv.org/abs/1609.07450) (2016). A new paper describing the assignment relaxation, the new local search and an extensive experimental study is in preparation.
* G. Farnebäck, [Go Game Graphs](https://github.com/GunnarFarneback/GoGameGraphs): large directed longest-path instances with published bounds.

## License and acknowledgements
Apache License 2.0 (see `LICENSE`). Research supported in part by PAPIIT IA106316, UNAM.
