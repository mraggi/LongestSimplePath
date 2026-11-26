# Methods, presets and options of `lsp solve`

## Methods (`--method`)

* `classic` — the 2018 schedule (paper v3 + code): 10 parameter sets; for each, forward DFS from 5 starts (ranked/random alternating) + erase 16 + backward DFS, then a PTO search; the same backward; a final PTO search. DFS runs stop after 0.1 s without improvement, PTO searches after 4 s without improvement (both configurable, also as work counts). It does **not** reproduce the 2018 reporting bugs (see `CHANGELOG.md`); everything else is faithful, including the quirks that can be switched off (`symmetric_in`, `in_score_typo`, `faithful_backward_greedy`, `largest_wcc`).
* `ils` — anytime iterated local search built from the same components: DFS construction, PTO improvement, then repeatedly keep a random window of the incumbent (`ils_keep` fraction), re-extend it by DFS both ways and re-run the PTO search; accept if not worse. Runs until `--time` (or `ils_rounds` rounds). With `ap_init=1` it first solves the assignment relaxation, patches its cycles into a path, starts from that path if it is better than the DFS start, and stops as soon as a path reaches the certified bound.
* `rdfs` — baseline: randomized DFS with random adjacency order and restarts from random vertices after `rdfs_steps` (default 10n) steps without improvement.
* `greedy` — baseline: repeated randomized Warnsdorff greedy (choose the free neighbour with fewest free continuations; dead ends last; heaviest arc first on weighted graphs), extending the tail and then the head.
* `gamm` — baseline: the mutation-based GA (GAMM) of Portugal, Antunes & Rocha (IEEE SMC 2010), re-implemented from the paper: degree-biased random initial paths grown in both directions, mutation = cut the path at a distance given by the adaptive perturbation pressure and regrow it (avoiding the old successor) at either end, duplicate removal, elitism (`ga_elite`), roulette selection, restart after `ga_stagnation` generations without improvement. Adaptation to digraphs: forward regrowth uses out-arcs/out-degrees, regrowth at the start uses in-arcs/in-degrees (the undirected original flips the path instead).

Rooted variant (`--source s`): paths must start at s. The graph is restricted to the vertices reachable from s, the in-arcs of s are removed, the DFS only starts at s, the DP only extends values reachable from s, and in the assignment relaxation the dummy vertex is forced to choose s.

## Presets

`--preset newap` (**the default**, with `--time 10`) = `--preset new` + `ap_init=1`. `--preset new` (the iterated local search alone): `method=ils`, `border=dfs`, `heavy_first=1`, `symmetric_in=1`, `in_score_typo=0`, `faithful_backward_greedy=0`, `largest_wcc=0`, `skip_block=1`, `min_free=1`, `dfs_steps=100000` without time rule, `pto_passes=3` without time rule. `--preset 2018`: method `classic` with the values of the table below, faithful to the 2018 version of this code (commit 7bc5804). A preset resets all algorithmic options but keeps `--time`, `--seed` and `--source`; put it before any `--set`.

## Options (`--set key=value`)

The column gives the value under `--preset 2018`; the presets `new` and `newap` change only the keys listed above.

| key | 2018 value | meaning |
|---|---|---|
| `largest_wcc` | 1 | keep only the largest weakly connected component (2018) |
| `rank` / `degree_rules` | 1 / 1 | rules 1 and 2–3 of the vertex orders |
| `score` | `low` | `low`, `high`, `random`, `none` |
| `symmetric_in` | 0 | 1 = mirror the out-order for the in-order (2018 sorts in-rank ascending) |
| `in_score_typo` | 1 | 1 = 2018 in-score (a3 added to a7) |
| `weighted_scores` | 0 | walk weights instead of walk counts |
| `heavy_first` | 0 | sort every adjacency list by arc weight (heaviest first), the vertex order breaking ties; no effect with unit weights |
| `param_restarts`, `params` | 10, `2018` | parameter sets: `2018` (P0, P1, random), `fixed`, `random` |
| `fwd_starts`, `bwd_starts`, `erase` | 5, 3, 16 | DFS starts and vertices erased before searching the other direction |
| `backward`, `other_direction`, `ranked_starts` | 1, 1, 1 | |
| `faithful_backward_greedy` | 1 | 2018: extend the back greedily before a backward DFS |
| `dfs_seconds`, `dfs_steps` | 0.1, 0 | DFS stops after this much time / these many steps without improvement |
| `pto`, `pto_each_dfs` | 1, 1 | PTO phase at all; after each DFS (else only at the end) |
| `pto_mode` | `open` | `open` (edge opening) or `reimpose` (fresh random PTO + imposition only, paper §3.2.1) |
| `pto_seconds`, `pto_passes` | 4, 0 | PTO search stops after this much time / opening passes without improvement |
| `pto_incremental` | 1 | incremental DP (0 = recompute from position 0; speed ablation) |
| `skip_block`, `min_free` | 4, 4 | blocks of size ≤ skip_block or with fewer than min_free free vertices are not opened |
| `border` | `first` | ordering of the free vertices B: `first` (2018 heuristic sort), `uniform`, `shuffle`, `dfs` (reverse post-order of a DFS of D[B]) |
| `sort_iters`, `sort_iters_per_free`, `sweep` | 12000, 0, 1 | heuristic-sort effort |
| `random_gaps`, `reshuffle_each_gap`, `try_reverse` | 1, 0, 0 | gap order; reshuffle B per gap; also try reverse(B) |
| `free_first` | 0 | emulate 2018 bug B1 (first path vertex treated as free in the opener; values may decrease; the best path seen is kept) |
| `sideways`, `sideways_max` | 0, 1000 | probability of accepting equal-weight openings with a different vertex set (plateau moves; tried, not helpful, off) |
| `ils_keep`, `ils_dfs_steps` | 0.5, 0 | ILS perturbation |
| `init_path` | (none) | ILS: file with a starting path (original vertex ids, whitespace-separated); used if it is a valid path (starting at `--source` if given) better than the DFS start |
| `ap_init` | 0 | ILS: start from the patched assignment relaxation (preset `newap`) |
| `ils_rounds` | 0 | ILS: stop after this many rounds (0 = only the time limit); for deterministic runs |
| `rdfs_steps` | 0 (=10n) | rdfs restart rule |
| `ga_population`, `ga_elite`, `ga_stagnation` | 100, 0.05, 15 | GAMM |

Determinism: with `--time 0` and only step-based stopping rules (`dfs_seconds=0 dfs_steps=K pto_seconds=0 pto_passes=K`) a run is a deterministic function of the graph, the options and `--seed` (unit test `WorkBasedRunsAreDeterministic`). Time-based rules make the amount of work machine dependent.

## Output of `lsp solve`

One JSON object: `value`, `num_vertices_path`, `check` (`ok` or the reason the internal checker rejected the path), `work_n`/`work_m` (graph actually searched), `read_seconds`, `preprocess_seconds`, `total_seconds`, counters (`dfs_steps`, `openings`, `open_passes`, `open_improvements`, `dp_vertices`, `sort_ops`, `pto_searches`, `restarts`), all `options`, `trajectory` = list of `[seconds since start of solve, value]` at every improvement, the assignment relaxation (`ap_solved`, `ap_certified`, `ap_bound_global`, `ap_bound` = upper bound on every path if certified and global, `ap_cycles`, `ap_patched`, `ap_seconds`, `patch_seconds`, ...), `proven_optimal`, and `path` with `--print-path` (vertex ids of the input file). Exit status 1 if the internal check fails.

`lsp ap GRAPH [--source s] [--path-out FILE]` solves only the relaxation and the patching and prints the bound and the patched path's weight.

`scripts/pathcheck.py` re-checks a path with a pure-Python checker sharing no code with the solver; `scripts/verify_ap_bound.py` recomputes the assignment bound with SciPy.
