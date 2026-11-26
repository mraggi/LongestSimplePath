// Strongly / weakly connected components and the ranks of the condensation.
#pragma once

#include "lsp/graph.hpp"

namespace lsp {

struct SccInfo
{
    // comp[v] = id of the SCC of v. Ids follow a topological order of the condensation
    // (every arc between different SCCs goes from a smaller to a larger id).
    std::vector<std::int32_t> comp;
    std::int32_t num_components = 0;
    // Block layout: members of component c are block[start[c] .. start[c+1]).
    std::vector<std::int32_t> start;
    std::vector<Vertex> block;
    // rank_out[c]: number of arcs of a longest path of the condensation starting at c.
    // rank_in[c]:  number of arcs of a longest path of the condensation ending at c.
    std::vector<std::int32_t> rank_out;
    std::vector<std::int32_t> rank_in;

    [[nodiscard]] std::int32_t size(std::int32_t c) const { return start[c + 1] - start[c]; }
    [[nodiscard]] std::int32_t rank_out_of(Vertex v) const { return rank_out[comp[v]]; }
    [[nodiscard]] std::int32_t rank_in_of(Vertex v) const { return rank_in[comp[v]]; }
};

SccInfo compute_scc(const Digraph& g); // iterative Tarjan, O(n + m)

// Weakly connected components: wcc[v] in [0, count).
std::vector<std::int32_t> weak_components(const Digraph& g, std::int32_t& count);

// Vertices reachable from s (including s).
std::vector<char> reachable_from(const Digraph& g, Vertex s);

} // namespace lsp
