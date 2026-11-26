#include "lsp/graph.hpp"

#include <stdexcept>

namespace lsp {

Digraph::Digraph(Vertex n, std::vector<Edge> edges) : n_(n)
{
    if (n < 0)
        throw std::invalid_argument("negative number of vertices");
    std::erase_if(edges, [](const Edge& e) { return e.from == e.to; });
    for (const auto& e : edges)
        if (e.from < 0 || e.from >= n || e.to < 0 || e.to >= n)
            throw std::invalid_argument("edge endpoint out of range");
    std::ranges::sort(edges, [](const Edge& a, const Edge& b) {
        if (a.from != b.from)
            return a.from < b.from;
        if (a.to != b.to)
            return a.to < b.to;
        return a.weight > b.weight; // keep the heaviest parallel arc first
    });
    // merge parallel arcs (keep max weight)
    std::vector<Edge> uniq;
    uniq.reserve(edges.size());
    for (const auto& e : edges)
        if (uniq.empty() || uniq.back().from != e.from || uniq.back().to != e.to)
            uniq.push_back(e);

    const auto un = static_cast<std::size_t>(n);
    out_off_.assign(un + 1, 0);
    in_off_.assign(un + 1, 0);
    for (const auto& e : uniq)
    {
        ++out_off_[static_cast<std::size_t>(e.from) + 1];
        ++in_off_[static_cast<std::size_t>(e.to) + 1];
        if (e.weight != 1)
            unit_weight_ = false;
    }
    for (std::size_t i = 0; i < un; ++i)
    {
        out_off_[i + 1] += out_off_[i];
        in_off_[i + 1] += in_off_[i];
    }
    out_arcs_.resize(uniq.size());
    in_arcs_.resize(uniq.size());
    std::vector<std::size_t> po(out_off_.begin(), out_off_.end() - 1);
    std::vector<std::size_t> pi(in_off_.begin(), in_off_.end() - 1);
    for (const auto& e : uniq) // sorted by (from,to): out-lists sorted by head
        out_arcs_[po[static_cast<std::size_t>(e.from)]++] = {e.to, e.weight};
    // in-lists sorted by tail: iterate by tail order
    for (const auto& e : uniq)
        in_arcs_[pi[static_cast<std::size_t>(e.to)]++] = {e.from, e.weight};
}

std::optional<Weight> Digraph::arc_weight(Vertex u, Vertex v) const
{
    if (u < 0 || u >= n_ || v < 0 || v >= n_)
        return std::nullopt;
    auto arcs = out(u);
    auto it = std::ranges::lower_bound(arcs, v, {}, &Arc::v);
    if (it != arcs.end() && it->v == v)
        return it->w;
    return std::nullopt;
}

std::vector<Edge> Digraph::edges() const
{
    std::vector<Edge> es;
    es.reserve(num_edges());
    for (Vertex u = 0; u < n_; ++u)
        for (const auto& a : out(u))
            es.push_back({u, a.v, a.w});
    return es;
}

Digraph Digraph::induced(const std::vector<char>& keep, std::vector<Vertex>& old_of_new) const
{
    std::vector<Vertex> new_of_old(static_cast<std::size_t>(n_), kNoVertex);
    old_of_new.clear();
    for (Vertex v = 0; v < n_; ++v)
        if (keep[static_cast<std::size_t>(v)])
        {
            new_of_old[static_cast<std::size_t>(v)] = static_cast<Vertex>(old_of_new.size());
            old_of_new.push_back(v);
        }
    std::vector<Edge> es;
    for (Vertex u = 0; u < n_; ++u)
    {
        const Vertex nu = new_of_old[static_cast<std::size_t>(u)];
        if (nu == kNoVertex)
            continue;
        for (const auto& a : out(u))
        {
            const Vertex nv = new_of_old[static_cast<std::size_t>(a.v)];
            if (nv != kNoVertex)
                es.push_back({nu, nv, a.w});
        }
    }
    Digraph h(static_cast<Vertex>(old_of_new.size()), std::move(es));
    h.labels.reserve(old_of_new.size());
    for (Vertex v : old_of_new)
        h.labels.push_back(label(v));
    return h;
}

Digraph Digraph::without_in_arcs(Vertex v) const
{
    std::vector<Edge> es;
    es.reserve(num_edges());
    for (Vertex u = 0; u < n_; ++u)
        for (const auto& a : out(u))
            if (a.v != v)
                es.push_back({u, a.v, a.w});
    Digraph h(n_, std::move(es));
    h.labels = labels;
    return h;
}

PathCheck check_path(const Digraph& g, const std::vector<Vertex>& vertices, std::optional<Weight> claimed,
                     std::optional<Vertex> required_source)
{
    PathCheck res;
    if (vertices.empty())
    {
        res.message = "empty path";
        return res;
    }
    std::vector<char> seen(static_cast<std::size_t>(g.num_vertices()), 0);
    for (std::size_t i = 0; i < vertices.size(); ++i)
    {
        const Vertex v = vertices[i];
        if (v < 0 || v >= g.num_vertices())
        {
            res.message = "vertex out of range at position " + std::to_string(i);
            return res;
        }
        if (seen[static_cast<std::size_t>(v)])
        {
            res.message = "repeated vertex " + std::to_string(v) + " at position " + std::to_string(i);
            return res;
        }
        seen[static_cast<std::size_t>(v)] = 1;
        if (i > 0)
        {
            auto w = g.arc_weight(vertices[i - 1], v);
            if (!w)
            {
                res.message = "missing arc " + std::to_string(vertices[i - 1]) + "->" + std::to_string(v);
                return res;
            }
            res.recomputed += *w;
        }
    }
    if (required_source && vertices.front() != *required_source)
    {
        res.message = "path does not start at the required source";
        return res;
    }
    if (claimed && *claimed != res.recomputed)
    {
        res.message = "claimed value " + std::to_string(*claimed) + " != recomputed " + std::to_string(res.recomputed);
        return res;
    }
    res.ok = true;
    return res;
}

} // namespace lsp
