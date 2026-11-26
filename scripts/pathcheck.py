"""Independent path checker (pure Python, shares no code with the C++ implementation).

check(graph_file, path, source=None) -> (ok, value, message)
Reads the 0-based edge-list format ("n m" header, "u v [w]" lines, '#' comments). Parallel
arcs keep the maximum weight, self-loops are ignored (same convention as the solver).
"""
import functools


@functools.lru_cache(maxsize=8)
def load(graph_file):
    w = {}
    n = None
    with open(graph_file) as f:
        for line in f:
            if not line.strip() or line[0] in "#%":
                continue
            p = line.split()
            if n is None:
                n = int(p[0])
                continue
            u, v = int(p[0]), int(p[1])
            x = int(p[2]) if len(p) > 2 else 1
            if u != v and ((u, v) not in w or x > w[(u, v)]):
                w[(u, v)] = x
    return n, w


def check(graph_file, path, source=None):
    n, w = load(graph_file)
    if not path:
        return False, 0, "empty path"
    if any(not (0 <= v < n) for v in path):
        return False, 0, "vertex out of range"
    if len(set(path)) != len(path):
        return False, 0, "repeated vertex"
    total = 0
    for a, b in zip(path, path[1:]):
        if (a, b) not in w:
            return False, 0, f"missing arc {a}->{b}"
        total += w[(a, b)]
    if source is not None and path[0] != source:
        return False, total, "wrong source"
    return True, total, "ok"


if __name__ == "__main__":
    import sys
    p = [int(x) for x in open(sys.argv[2]).read().split()]
    print(check(sys.argv[1], p, int(sys.argv[3]) if len(sys.argv) > 3 else None))
