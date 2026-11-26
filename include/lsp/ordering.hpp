// Vertex scores and the vertex orders that drive the depth-first search (paper §3.1).
#pragma once

#include "lsp/structure.hpp"

#include <array>

namespace lsp {

using Params = std::array<double, 8>; // (a1..a4) out-score, (a5..a8) in-score

enum class ScoreRule
{
    Low,    // smaller score first (2018 default: "save the good vertices for later")
    High,   // larger score first
    Random, // uniformly random key, redrawn for every parameter set
    None,   // ignore the score (ties broken by vertex index)
};

struct OrderingOptions
{
    bool use_rank = true;          // rule 1: larger out-rank first
    bool use_degree_rules = true;  // rules 2-3: sinks last, in-degree-1 vertices first
    ScoreRule score = ScoreRule::Low;
    bool symmetric_in = false;     // false: 2018 (ascending in-rank in the in-order); true: mirror image
    bool in_score_typo = true;     // true: 2018 in-score (a3 added to a7)
    bool weighted_scores = false;  // use total walk weight instead of walk counts
    bool heavy_arcs_first = false; // sort each adjacency list by arc weight (descending) first
};

struct Scores
{
    std::vector<double> out;
    std::vector<double> in;
};
Scores compute_scores(const Digraph& g, const Params& a, const OrderingOptions& opt);

// Adjacency lists re-sorted by a ranking of the vertices.
struct OrderedAdjacency
{
    std::vector<std::size_t> off;
    std::vector<Arc> arcs;
    [[nodiscard]] std::span<const Arc> operator[](Vertex v) const
    {
        return {arcs.data() + off[static_cast<std::size_t>(v)], arcs.data() + off[static_cast<std::size_t>(v) + 1]};
    }
};

struct SearchOrders
{
    std::vector<Vertex> out_order; // all vertices sorted by the out-order (forward DFS starts)
    std::vector<Vertex> in_order;  // all vertices sorted by the in-order (backward DFS starts)
    OrderedAdjacency out_adj;      // out-lists sorted by the out-order
    OrderedAdjacency in_adj;       // in-lists sorted by the in-order
};

SearchOrders make_orders(const Digraph& g, const SccInfo& scc, const Params& a, const OrderingOptions& opt, Rng& rng);

// Adjacency in the input order (used by baselines).
OrderedAdjacency plain_out_adjacency(const Digraph& g);
OrderedAdjacency plain_in_adjacency(const Digraph& g);

// The parameter sets used by code_2018 (GetParams): P0, P1 ("trained"), then random integers in [-10, 89].
Params params_2018(int i, Rng& rng);

} // namespace lsp
