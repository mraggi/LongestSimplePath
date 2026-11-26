// Weighted digraph in compressed sparse row form (out- and in-adjacency).
#pragma once

#include "lsp/core.hpp"

#include <optional>
#include <string>
#include <vector>

namespace lsp {

struct Edge
{
    Vertex from;
    Vertex to;
    Weight weight;
};

struct Arc
{
    Vertex v;   // the other endpoint (head for out-arcs, tail for in-arcs)
    Weight w;
};

class Digraph
{
public:
    Digraph() = default;

    // Builds the digraph. Self-loops are dropped (a simple path never uses them).
    // Parallel arcs u->v are merged keeping the maximum weight.
    Digraph(Vertex n, std::vector<Edge> edges);

    [[nodiscard]] Vertex num_vertices() const { return n_; }
    [[nodiscard]] std::size_t num_edges() const { return out_arcs_.size(); }

    [[nodiscard]] std::span<const Arc> out(Vertex u) const
    {
        return {out_arcs_.data() + out_off_[u], out_arcs_.data() + out_off_[u + 1]};
    }
    [[nodiscard]] std::span<const Arc> in(Vertex v) const
    {
        return {in_arcs_.data() + in_off_[v], in_arcs_.data() + in_off_[v + 1]};
    }
    [[nodiscard]] std::size_t outdegree(Vertex u) const { return out_off_[u + 1] - out_off_[u]; }
    [[nodiscard]] std::size_t indegree(Vertex v) const { return in_off_[v + 1] - in_off_[v]; }

    // Weight of arc u->v if present (binary search; out-lists are sorted by head).
    [[nodiscard]] std::optional<Weight> arc_weight(Vertex u, Vertex v) const;

    [[nodiscard]] bool is_unit_weight() const { return unit_weight_; }
    [[nodiscard]] std::vector<Edge> edges() const;

    // Optional external labels (e.g. the original vertex ids of a SNAP file).
    std::vector<std::string> labels;
    [[nodiscard]] std::string label(Vertex v) const
    {
        return labels.empty() ? std::to_string(v) : labels[static_cast<std::size_t>(v)];
    }

    // Induced subgraph on the vertices with keep[v] = true. Returns the new graph and
    // fills old_of_new with the original index of every new vertex. Labels are carried.
    [[nodiscard]] Digraph induced(const std::vector<char>& keep, std::vector<Vertex>& old_of_new) const;

    // Same digraph without the arcs entering `v` (used for the rooted variant).
    [[nodiscard]] Digraph without_in_arcs(Vertex v) const;

private:
    Vertex n_ = 0;
    std::vector<std::size_t> out_off_{0};
    std::vector<std::size_t> in_off_{0};
    std::vector<Arc> out_arcs_;
    std::vector<Arc> in_arcs_;
    bool unit_weight_ = true;
};

// A path as a list of vertices plus its weight.
struct Path
{
    std::vector<Vertex> vertices;
    Weight value = 0;

    [[nodiscard]] std::size_t num_edges() const { return vertices.empty() ? 0 : vertices.size() - 1; }
};

// Independent verification: simple, every consecutive pair is an arc, weight matches.
struct PathCheck
{
    bool ok = false;
    Weight recomputed = 0;
    std::string message;
};
PathCheck check_path(const Digraph& g, const std::vector<Vertex>& vertices, std::optional<Weight> claimed = std::nullopt,
                     std::optional<Vertex> required_source = std::nullopt);

} // namespace lsp
