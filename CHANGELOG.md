# Changelog

## 2026: C++20 rewrite

The code was rewritten from scratch in C++20. The 2016–2018 version is in the git history (last commit 7bc5804).

### New
* **Assignment relaxation.** Gives an upper bound valid for arbitrary weights and a starting path, obtained by patching the cycles of the assignment. It is solved by an ε-scaling auction and its optimality is checked with an exact integer dual certificate. The solver stops with a proof of optimality when the path reaches the bound (`--preset newap`).
* **Pseudo-topological local search.** The free vertices of a strongly connected component are ordered by the reverse post-order of a depth-first search started next to the gap, instead of the randomized "heuristic sort". This is the main improvement on geometric graphs, grids and Go game graphs.
* **Anytime iterated local search** (`--preset new`):
  * deterministic, work-based stopping rules;
  * runs until a time limit;
  * heavier arcs first on weighted graphs;
  * every non-trivial component is opened;
  * the whole graph is searched, not only its largest weakly connected component.
* **Rooted variant** (`--source s`).
* **Command-line tools** `lsp` and `lsp_check`, reading edge lists, Go Game Graphs adjacency lists, SNAP and DIMACS files.
* **Verification and baselines:**
  * an independent path checker (`lsp_check`, `scripts/pathcheck.py`);
  * an independent SciPy check of the bound (`scripts/verify_ap_bound.py`);
  * unit tests;
  * baselines: randomized DFS, Warnsdorff greedy, the mutation GA of Portugal, Antunes & Rocha (2010);
  * a CP-SAT hybrid (`scripts/hybrid.py`).
* Every component can be switched off for experiments (`--set key=value`, see `docs/OPTIONS.md`).

### Bugs of the 2018 code fixed
1. The first vertex of every path produced by the improvement phase was dropped from the returned path, although its arc was counted in the reported weight. The reported value was therefore too high by one arc.
2. A stale parent pointer in the dynamic program could produce invalid paths (consecutive vertices not joined by an arc). This happened in all runs on a graph with negative weights, and on some positively weighted random graphs.
3. The incrementally maintained maximum ignored the unchanged prefix of the order in one case.

Every returned path is now verified (`check` field of the output).

## 2015–2018
Original version, written for the Oracle MDC 2015 coding challenge and described in arXiv:1609.07450.
