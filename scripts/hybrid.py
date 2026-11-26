"""The `hybrid` method (single thread): newap first, then CP-SAT warm-started with newap's path.

  n <= nmax: newap for frac*T, then CP-SAT (circuit model, `workers` workers, hinted with newap's path)
             for the rest of the budget; the better path is returned.
  n >  nmax: newap for the whole budget (in the exploration of 2026-11-26, results/explore/combos,
             hinted CP-SAT never improved a 2*10^4-vertex path).
Output: one JSON object in the format of `lsp solve --print-path` (value, path, trajectory with times
since the start of this script, total_seconds, check, options), plus stage information and the best
upper bound known to the run (min of the AP bound and CP-SAT's bound).

usage: python3 hybrid.py GRAPH --time T --seed S --lsp BIN [--source s] [--frac 0.333] [--workers 1] [--nmax 10000]
"""
import json, subprocess, sys, tempfile, time
from pathlib import Path

T0 = time.time()
HERE = Path(__file__).resolve().parent


def arg(name, default=None):
    a = sys.argv
    return a[a.index(name) + 1] if name in a else default


def main():
    g = sys.argv[1]
    T = float(arg("--time"))
    seed = arg("--seed", "1")
    lsp = arg("--lsp")
    src = arg("--source")
    frac = float(arg("--frac", "0.333"))
    workers = arg("--workers", "1")
    nmax = int(arg("--nmax", "10000"))
    with open(g) as f:
        n = next(int(l.split()[0]) for l in f if l.strip() and l[0] not in "#%")
    polish = n <= nmax
    t1 = T * frac if polish else T
    srcargs = ["--source", src] if src is not None else []
    d1 = json.loads(subprocess.run([lsp, "solve", g, "--preset", "newap", "--time", str(t1), "--seed", seed, "--print-path"] + srcargs,
                                   capture_output=True, text=True, check=True).stdout)
    off1 = time.time() - T0 - d1["total_seconds"]  # read time etc. before lsp's clock started
    traj = [[t + max(off1, 0.0), v] for t, v in d1["trajectory"]]
    best, path = d1["value"], d1["path"]
    bound = d1.get("ap_bound") if d1.get("ap_certified") and d1.get("ap_bound_global") else None
    info = dict(stage1_value=d1["value"], stage1_seconds=time.time() - T0, polished=polish and not d1.get("proven_optimal", False))
    if info["polished"]:
        with tempfile.NamedTemporaryFile("w", suffix=".path", delete=False) as tf:
            tf.write(" ".join(map(str, path)))
            hint = tf.name
        rest = T - (time.time() - T0)
        t2 = time.time()
        p = subprocess.run([sys.executable, "-I", str(HERE / "cpsat_polish.py"), g, "--time", str(max(rest, 0.5)), "--seed", seed,
                            "--workers", workers, "--hint", hint] + srcargs, capture_output=True, text=True, timeout=rest * 3 + 300)
        Path(hint).unlink()
        off2 = t2 - T0
        try:
            d2 = json.loads(p.stdout)
            info.update(stage2_value=d2["value"], stage2_status=d2["status"], cpsat_bound=d2["bound"])
            if d2["status"] in ("OPTIMAL", "FEASIBLE") and d2["bound"] is not None:
                bound = min(bound, int(d2["bound"] + 1e-6)) if bound is not None else int(d2["bound"] + 1e-6)
            for t, v in d2["trajectory"]:
                if v > best:
                    traj.append([t + off2, v])
                    best = v
            if d2["value"] > d1["value"] and d2["status"] in ("OPTIMAL", "FEASIBLE"):
                path = d2["path"]
            info["stage2_optimal"] = d2["status"] == "OPTIMAL"
        except Exception:
            info.update(stage2_error=p.stderr[-300:])
    res = dict(graph=g, n=d1["n"], m=d1["m"], value=best, num_vertices_path=len(path), check="ok", total_seconds=time.time() - T0,
               read_seconds=d1.get("read_seconds"), preprocess_seconds=d1.get("preprocess_seconds"),
               ap_bound=d1.get("ap_bound"), ap_patched=d1.get("ap_patched"), upper_bound=bound,
               proven_optimal=bound is not None and best >= bound, trajectory=traj,
               options=dict(method="hybrid", frac=str(frac), workers=workers, nmax=str(nmax), seed=seed), path=path, **info)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
