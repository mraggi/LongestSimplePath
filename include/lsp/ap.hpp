// Assignment relaxation (AP) of the longest simple path problem, and Karp-style patching of the
// AP cycles into one path.
//
// The relaxation. Add a dummy vertex D. Every vertex v chooses a successor among its out-neighbours,
// itself (v is skipped) and D (v ends the path); D chooses the first vertex of the path, or itself
// (empty path). The choices must form a permutation of V+D. Its cycle through D is a simple path,
// and its other non-trivial cycles are vertex-disjoint cycles of the digraph. The weight of the
// permutation is the total weight of the digraph arcs it uses. Every simple path P is the D-cycle of
// such a permutation (the other vertices skipped), so the maximum weight of a permutation is an
// upper bound on the weight of every simple path, for arbitrary arc weights. With a source s: D->s is
// forced, arcs into s and the skip of s are not allowed, and the bound holds for paths starting at s.
// This is the LP relaxation (without subtour constraints) of the circuit model that CP-SAT solves.
#pragma once

#include "lsp/graph.hpp"

namespace lsp {

struct ApSolution
{
    bool solved = false;                      // false if the deadline expired first
    bool certified = false;                   // optimality certified by the auction prices (dual bound)
    Weight bound = 0;                         // optimal AP value: upper bound on every simple path
    std::vector<Vertex> path;                 // the cycle through D, without D
    std::vector<std::vector<Vertex>> cycles;  // the other cycles of length >= 2, each in cycle order
    std::int64_t bids = 0;                    // auction work
    int phases = 0;                           // epsilon-scaling phases
};

// Solves the AP exactly with an epsilon-scaling forward auction (Bertsekas) on integer benefits
// scaled by (#rows + 1), so that the final phase (epsilon = 1) is optimal. The optimality is also
// checked independently: the prices give a dual (upper) bound, and `certified` is set when it
// differs from the assignment value by less than one unscaled unit.
ApSolution solve_ap(const Digraph& g, Vertex source = kNoVertex, const Deadline& deadline = Deadline::never());

struct PatchStats
{
    std::int32_t merged = 0;      // cycles merged with a gain (append, prepend, 2-exchange insertion)
    std::int32_t lossy = 0;       // cycles merged by cutting off a prefix or suffix of the path
    std::int32_t left = 0;        // cycles not merged
    std::int32_t extended = 0;    // vertices added by the final greedy extension
};

// Turns an AP solution into one simple path of g (starting at `source` if given):
// (1) merge cycles into the path while that gains: append at the end, prepend at the start, or
//     insert between consecutive path vertices a, a' (arcs a->c', c->a', where c->c' is a cycle arc);
// (2) merge a cycle by cutting off a suffix or a prefix of the path, if the cycle is worth more
//     (best merge over all cycles; then (1) again);
// (3) extend both ends greedily (heaviest arc, then largest degree).
// Prepending and prefix cuts are disabled when the path must start at `source`.
std::vector<Vertex> patch_ap(const Digraph& g, const ApSolution& ap, Vertex source = kNoVertex, PatchStats* stats = nullptr);

} // namespace lsp
