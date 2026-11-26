// Minimal use of the library: build a weighted digraph and search for a heavy simple path.
#include "lsp/solver.hpp"

#include <iostream>

int main()
{
    // vertices 0..5; arcs (from, to, weight); weights may be negative
    lsp::Digraph g(6, {{0, 1, 3}, {1, 2, 2}, {2, 0, 4}, {2, 3, 1}, {3, 4, 5}, {1, 4, 1}, {4, 5, 2}, {5, 3, -1}});

    lsp::SolverOptions opt;
    opt.apply_preset("newap"); // assignment opening + pseudo-topological local search
    opt.time_limit = 1;        // seconds (it also stops as soon as the path is proven optimal)
    opt.seed = 1;
    // opt.source = 0;         // rooted variant: paths must start at vertex 0

    const lsp::SolveResult res = lsp::solve(g, opt);
    std::cout << "weight " << res.best.value << (res.proven_optimal ? " (proven optimal)" : "") << "\npath";
    for (lsp::Vertex v : res.best.vertices)
        std::cout << ' ' << v;
    std::cout << '\n';
    if (res.ap.certified && res.ap.bound_global)
        std::cout << "upper bound (assignment relaxation): " << res.ap.bound << '\n';
}
