"""Independent check of the assignment-relaxation (AP) bound reported by `lsp` (code/src/ap.cpp).

Solves the same relaxation with scipy's sparse LAPJVsp (min_weight_full_bipartite_matching), which shares
no code with the C++ auction, and prints the optimum. The relaxation (paper, Section 3.3): rows and columns
are the vertices plus a dummy D; each vertex v picks one successor among its out-neighbours (value w(vu)),
itself (v skipped, value 0) or D (the path ends at v, value 0); D picks the first vertex of the path.
All arcs are kept, whatever their sign. Rooted variant (--source s): only vertices reachable from s,
D -> s forced, no arc enters s, s is not skipped.

usage: python3 verify_ap_bound.py GRAPH [--source s]
"""
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathcheck import load  # noqa: E402


def main():
    a = sys.argv
    g = a[1]
    src = int(a[a.index("--source") + 1]) if "--source" in a else None
    t0 = time.time()
    n, w = load(g)
    keep = np.ones(n, bool)
    if src is not None:
        out = [[] for _ in range(n)]
        for (u, v) in w:
            out[u].append(v)
        keep[:] = False
        keep[src] = True
        q = deque([src])
        while q:
            u = q.popleft()
            for v in out[u]:
                if not keep[v]:
                    keep[v] = True
                    q.append(v)
        del out
    idx = -np.ones(n, np.int64)
    verts = np.flatnonzero(keep)
    idx[verts] = np.arange(len(verts))
    k = len(verts)
    D = k
    arcs = [(idx[u], idx[v], x) for (u, v), x in w.items() if keep[u] and keep[v] and v != src]
    W = max([abs(x) for _, _, x in arcs] + [1])
    big = W + 1  # cost of a choice of value x is big - x > 0; every row picks exactly one entry
    rows = [u for u, _, _ in arcs]
    cols = [v for _, v, _ in arcs]
    cost = [big - x for _, _, x in arcs]
    s = idx[src] if src is not None else None
    for v in range(k):
        rows.append(v); cols.append(D); cost.append(big)          # v ends the path
        if v != s:
            rows.append(v); cols.append(v); cost.append(big)      # v skipped
    if s is None:
        for v in range(k):
            rows.append(D); cols.append(v); cost.append(big)      # the path starts at v
        rows.append(D); cols.append(D); cost.append(big)          # empty path
    else:
        rows.append(D); cols.append(s); cost.append(big)
    M = csr_matrix((np.array(cost, float), (np.array(rows), np.array(cols))), shape=(k + 1, k + 1))
    t1 = time.time()
    ri, ci = min_weight_full_bipartite_matching(M)
    total = M[ri, ci].sum()
    value = (k + 1) * big - total
    print(f"n={n} kept={k} arcs={len(arcs)} AP optimum (scipy) = {int(round(value))}  (build {t1 - t0:.1f}s, solve {time.time() - t1:.1f}s)")


if __name__ == "__main__":
    main()
