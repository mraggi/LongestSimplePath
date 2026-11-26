// Depth-first enumeration of simple paths with a single in-place path (paper §3.1, Algorithm 2).
//
// The path is stored as a sequence seq[0..k) grown and shrunk at its end. For a forward search
// the sequence is the path itself and `adj` is the out-adjacency; for a backward search the
// sequence is the reversed path and `adj` is the in-adjacency (arc weights are the same).
// seq[0] is the root and is never removed.
#pragma once

#include "lsp/ordering.hpp"

#include <functional>

namespace lsp {

struct StopRule
{
    double seconds_without_improvement = 0; // <= 0: not used
    std::int64_t steps_without_improvement = 0; // <= 0: not used
};

struct DfsStats
{
    std::int64_t steps = 0;
    std::int64_t improvements = 0;
};

using ValueSink = std::function<void(Weight)>; // called on every improvement of the search's best

class InPlaceDfs
{
public:
    InPlaceDfs(Vertex n, const OrderedAdjacency& adj);

    // Sets the current sequence (must be a valid walk in `adj` without repeated vertices).
    void reset(const std::vector<Vertex>& seq);
    // Appends the first unused neighbour of the last vertex while possible.
    void greedy_extend();
    // One step of the enumeration (NextPath + greedy extension). False when finished.
    bool step();

    [[nodiscard]] Weight value() const { return value_; }
    [[nodiscard]] std::size_t size() const { return seq_.size(); }
    [[nodiscard]] std::vector<Vertex> sequence() const;

    // Runs the enumeration from the current state until the stop rule or the deadline fires,
    // or the enumeration is exhausted. Returns the best sequence seen (including the start).
    std::vector<Vertex> run(const StopRule& stop, const Deadline& deadline, DfsStats& stats, const ValueSink& sink);

private:
    struct Entry
    {
        Vertex v;
        std::int32_t idx; // index of v in adj[pred]
        Weight w;         // weight of the arc pred -> v
    };
    std::int32_t first_free(Vertex u, std::int32_t from) const;
    void push(Vertex u, std::int32_t i);
    void pop();
    void sync_best();

    const OrderedAdjacency& adj_;
    std::vector<Entry> seq_;
    std::vector<char> used_;
    Weight value_ = 0;

    // lazily synchronised copy of the best sequence
    std::vector<Entry> best_;
    Weight best_value_ = 0;
    bool best_dirty_ = false;
    std::size_t common_ = 0; // seq_[0..common_) == best_[0..common_) when !best_dirty_
};

// Convenience: reverse a vertex sequence.
std::vector<Vertex> reversed(std::vector<Vertex> v);

} // namespace lsp
