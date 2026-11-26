"""CP-SAT (circuit model) warm-started with a path (--hint FILE), used as the polishing stage of the
`hybrid` method (results/scripts/hybrid.py). Copied from results/explore/ap_patching/scripts/cpsat_hint.py.
Original docstring follows.

CP-SAT used as an anytime heuristic (same circuit model as exact_cpsat.py), with output in the
JSON format of `lsp solve` so that the experiment runner can treat it as one more method.

usage: python3 cpsat_heuristic.py GRAPH --time T --seed S --workers W [--source s]
The trajectory records (seconds since the start of this process, objective) at every improving
solution, so model construction time is included in the budget; --time limits the whole run
(the solver gets the time that remains after building the model).
"""
import json, sys, time

T0 = time.time()
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathcheck import load  # noqa: E402


def main():
    a = sys.argv
    g = a[1]
    tl = float(a[a.index("--time") + 1])
    seed = int(a[a.index("--seed") + 1]) if "--seed" in a else 1
    workers = int(a[a.index("--workers") + 1]) if "--workers" in a else 1
    src = int(a[a.index("--source") + 1]) if "--source" in a else None
    from ortools.sat.python import cp_model
    n, w = load(g)
    out = [[] for _ in range(n)]
    for (u, v), x in w.items():
        out[u].append(v)
    keep = [True] * n
    if src is not None:
        keep = [False] * n
        keep[src] = True
        st = [src]
        while st:
            u = st.pop()
            for v in out[u]:
                if not keep[v]:
                    keep[v] = True
                    st.append(v)
    m = cp_model.CpModel()
    D = n
    arcs, obj, lit = [], [], {}
    for (u, v), x in w.items():
        if not (keep[u] and keep[v]) or (src is not None and v == src):
            continue
        b = m.NewBoolVar("")
        arcs.append((u, v, b))
        lit[(u, v)] = b
        if x:
            obj.append(x * b)
    for v in range(n):
        if not keep[v]:
            continue
        loop = m.NewBoolVar("")
        arcs.append((v, v, loop))
        if src is None or v == src:
            arcs.append((D, v, m.NewBoolVar("")))
        arcs.append((v, D, m.NewBoolVar("")))
        if src is not None and v == src:
            m.Add(loop == 0)
    m.AddCircuit(arcs)
    if "--hint" in a:
        hp = [int(x) for x in open(a[a.index("--hint") + 1]).read().split()]
        hs = {hp[i]: hp[i + 1] for i in range(len(hp) - 1)}
        inp = set(hp)
        for (u, v, b) in arcs:
            if u == D:
                val = v == hp[0]
            elif v == D:
                val = u == hp[-1]
            elif u == v:
                val = u not in inp
            else:
                val = hs.get(u) == v
            m.AddHint(b, int(val))
    m.Maximize(sum(obj))
    build = time.time() - T0
    traj = []

    class CB(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()

        def on_solution_callback(self):
            v = int(round(self.ObjectiveValue()))
            if not traj or v > traj[-1][1]:
                traj.append([time.time() - T0, v])

    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = max(0.1, tl - (time.time() - T0))
    s.parameters.num_workers = workers
    s.parameters.random_seed = seed
    status = s.Solve(m, CB())
    path, value = [], 0
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        succ = {u: v for (u, v), b in lit.items() if s.Value(b)}
        heads = set(succ.values())
        start = [u for u in succ if u not in heads]
        if start:
            u = start[0]
            path = [u]
            while u in succ:
                u = succ[u]
                path.append(u)
        elif src is not None:
            path = [src]
        value = int(round(s.ObjectiveValue()))
    if not path:  # no solution found: report the trivial path
        path = [src if src is not None else 0]
        value = 0
        if not traj:
            traj = [[time.time() - T0, 0]]
    res = dict(graph=g, n=n, m=len(w), value=value, num_vertices_path=len(path), check="ok", total_seconds=time.time() - T0,
               preprocess_seconds=build, status=s.StatusName(status), bound=s.BestObjectiveBound(), trajectory=traj,
               options=dict(method="cpsat", workers=str(workers), seed=str(seed)), path=path)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
