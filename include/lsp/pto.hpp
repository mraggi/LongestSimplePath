// Pseudo-topological orders, the longest-path DP on the forward subgraph, imposition of a path
// and the edge-opening move (paper §3.2).
#pragma once

#include "lsp/ordering.hpp"

namespace lsp {

enum class BOrderRule
{
    HeuristicFirst,   // 2018: sweep + random swaps with the FIRST earlier out-neighbour
    HeuristicUniform, // same, but the earlier out-neighbour is chosen uniformly (provably converges)
    ShuffleOnly,      // uniformly random order of B (no heuristic sort)
    DfsPostorder,     // new: reverse postorder of a randomized DFS of D[B], rooted first at the
                      // out-neighbours of the vertex before the gap (a weak PTO of D[B] in O(|B|+m_B))
};

struct OpenerOptions
{
    std::int32_t skip_block_size = 4; // blocks (SCCs) with at most this many vertices are not opened
    std::int32_t min_free = 4;        // blocks with fewer free (non-path) vertices are not opened
    BOrderRule border = BOrderRule::HeuristicFirst;
    std::int64_t sort_iters = 12000;  // random heuristic-sort operations per opening (2018: 12000)
    double sort_iters_per_free = 0;   // if > 0: use ceil(this * |B|) operations instead
    bool sweep = true;                // deterministic sweep before the random operations (2018)
    bool random_gap_order = true;     // gaps of a block are visited in random order (2018)
    bool reshuffle_each_gap = false;  // reshuffle B before every gap (2018: once per block)
    bool try_reverse = false;         // after each opening also try B reversed (strength-2 covering)
    bool free_first = false;          // emulate 2018 bug B1: the first path vertex is treated as free,
                                      // so it may be moved into B and the DP value may decrease
    double sideways = 0;              // probability of accepting an opening whose DP path has the same
                                      // weight but a different vertex set (plateau move)
};

struct OpenStats
{
    std::int64_t openings = 0;
    std::int64_t passes = 0;
    std::int64_t improvements = 0;
    std::int64_t dp_vertices = 0; // vertices re-evaluated by the incremental DP
    std::int64_t sort_ops = 0;
    std::int64_t sideways_moves = 0;
    bool last_pass_sideways = false;   // set by open_pass when it ended with a plateau move
};

class PseudoTopoOrder
{
public:
    // in_adj: in-lists used by the DP (order only affects ties). out_adj: out-lists used by the
    // heuristic sort ("first out-neighbour"). source != kNoVertex: rooted variant (paths must
    // start at source, which must have no in-arcs).
    PseudoTopoOrder(const Digraph& g, const SccInfo& scc, const OrderedAdjacency& in_adj,
                    const OrderedAdjacency& out_adj, Vertex source = kNoVertex);

    // Blocks in the SCC order of `scc`, uniformly random order inside every block.
    void randomize(Rng& rng);
    // Imposition T_P: the positions occupied by P are refilled with P in path order.
    void impose(const std::vector<Vertex>& path);

    // Optimum of the DP on the forward subgraph (evaluated lazily and incrementally).
    Weight value();
    Path best_path();
    [[nodiscard]] Vertex at(std::int32_t i) const { return order_[static_cast<std::size_t>(i)]; }
    [[nodiscard]] std::int32_t position(Vertex v) const { return pos_[static_cast<std::size_t>(v)]; }
    [[nodiscard]] const std::vector<Vertex>& order() const { return order_; }

    // One call of the 2018 eXtreme_edge_opener: tries the gaps of every eligible block and
    // returns true at the first strict improvement of the DP value.
    bool open_pass(const OpenerOptions& opt, Rng& rng, const Deadline& deadline, OpenStats& st);
    // Whether any block is eligible for opening (otherwise open_pass can never improve).
    [[nodiscard]] bool has_openable_block(const OpenerOptions& opt) const;

    bool incremental = true; // false: always re-evaluate the DP from position 0 (for ablation)
    void swap_positions(std::int32_t i, std::int32_t j);
    // Checks that every SCC occupies its block (strong PTO with the fixed block order).
    [[nodiscard]] bool is_strong_pto() const;
    // Reorders positions [b,c) by the reverse post-order of a randomized DFS of the induced
    // subgraph on the vertices there (public for testing).
    // `anchor`: the path vertex just before the gap (its free out-neighbours are the first roots).
    void order_free_dfs(std::int32_t b, std::int32_t c, Vertex anchor, Rng& rng);

private:
    void touch(std::int32_t i) { dirty_ = std::min(dirty_, i); }
    void evaluate(OpenStats* st = nullptr);
    std::vector<std::int32_t> path_positions();
    void transfer(std::int32_t a, std::int32_t b, std::int32_t c, std::int32_t h);
    void order_free(std::int32_t b, std::int32_t c, Vertex anchor, const OpenerOptions& opt, Rng& rng, OpenStats& st);
    std::int32_t earlier_out_neighbour(std::int32_t lo, std::int32_t r, BOrderRule rule, Rng& rng);

    const Digraph& g_;
    const SccInfo& scc_;
    const OrderedAdjacency& in_adj_;
    const OrderedAdjacency& out_adj_;
    Vertex source_;
    std::int32_t n_;
    std::vector<Vertex> order_;
    std::vector<std::int32_t> pos_;
    std::vector<Weight> dp_;
    std::vector<std::int32_t> parent_;
    std::vector<std::int32_t> prefbest_;
    std::int32_t dirty_ = 0;
    std::vector<std::int32_t> scratch_;
    std::vector<std::uint32_t> stamp_;
    std::uint32_t cur_stamp_ = 0;
    std::vector<Vertex> roots_;
    std::vector<Vertex> post_;
    struct DfsFrame
    {
        Vertex v;
        std::uint32_t k;
        std::uint32_t off;
    };
    std::vector<DfsFrame> dfs_stack_;
};

} // namespace lsp
