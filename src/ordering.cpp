#include "lsp/ordering.hpp"

#include <numeric>
#include <tuple>

namespace lsp {

namespace {

// W_k(v) = number (or total weight) of walks with k arcs starting at v, k = 1..4, following `adj`.
template <class AdjFn>
std::array<std::vector<double>, 4> walk_counts(Vertex n, AdjFn adj, bool weighted)
{
    const auto un = static_cast<std::size_t>(n);
    std::array<std::vector<double>, 4> cnt;
    std::array<std::vector<double>, 4> wt;
    for (auto& c : cnt)
        c.assign(un, 0.0);
    for (auto& c : wt)
        c.assign(un, 0.0);
    for (Vertex v = 0; v < n; ++v)
        for (const auto& a : adj(v))
        {
            cnt[0][static_cast<std::size_t>(v)] += 1.0;
            wt[0][static_cast<std::size_t>(v)] += static_cast<double>(a.w);
        }
    for (std::size_t k = 1; k < 4; ++k)
        for (Vertex v = 0; v < n; ++v)
        {
            double c = 0;
            double w = 0;
            for (const auto& a : adj(v))
            {
                const auto x = static_cast<std::size_t>(a.v);
                c += cnt[k - 1][x];
                w += static_cast<double>(a.w) * cnt[k - 1][x] + wt[k - 1][x];
            }
            cnt[k][static_cast<std::size_t>(v)] = c;
            wt[k][static_cast<std::size_t>(v)] = w;
        }
    return weighted ? wt : cnt;
}

} // namespace

Scores compute_scores(const Digraph& g, const Params& a, const OrderingOptions& opt)
{
    const Vertex n = g.num_vertices();
    const auto W = walk_counts(n, [&](Vertex v) { return g.out(v); }, opt.weighted_scores);
    const auto V = walk_counts(n, [&](Vertex v) { return g.in(v); }, opt.weighted_scores);
    Scores s;
    s.out.resize(static_cast<std::size_t>(n));
    s.in.resize(static_cast<std::size_t>(n));
    const double a3_in = opt.in_score_typo ? a[2] + a[6] : a[6];
    for (std::size_t v = 0; v < static_cast<std::size_t>(n); ++v)
    {
        s.out[v] = a[0] * W[0][v] + a[1] * W[1][v] + a[2] * W[2][v] + a[3] * W[3][v];
        s.in[v] = a[4] * V[0][v] + a[5] * V[1][v] + a3_in * V[2][v] + a[7] * V[3][v];
    }
    return s;
}

namespace {

using Key = std::tuple<long long, int, int, double, Vertex>;

OrderedAdjacency sorted_adjacency(Vertex n, const std::vector<std::int32_t>& rank_pos, bool use_out, const Digraph& g,
                                  bool heavy_first = false)
{
    OrderedAdjacency r;
    const auto un = static_cast<std::size_t>(n);
    r.off.assign(un + 1, 0);
    for (Vertex v = 0; v < n; ++v)
        r.off[static_cast<std::size_t>(v) + 1] = r.off[static_cast<std::size_t>(v)] + (use_out ? g.outdegree(v) : g.indegree(v));
    r.arcs.resize(r.off.back());
    for (Vertex v = 0; v < n; ++v)
    {
        auto src = use_out ? g.out(v) : g.in(v);
        auto first = r.arcs.begin() + static_cast<std::ptrdiff_t>(r.off[static_cast<std::size_t>(v)]);
        std::ranges::copy(src, first);
        if (!rank_pos.empty())
            std::sort(first, first + static_cast<std::ptrdiff_t>(src.size()), [&](const Arc& x, const Arc& y) {
                if (heavy_first && x.w != y.w)
                    return x.w > y.w;
                return rank_pos[static_cast<std::size_t>(x.v)] < rank_pos[static_cast<std::size_t>(y.v)];
            });
    }
    return r;
}

} // namespace

SearchOrders make_orders(const Digraph& g, const SccInfo& scc, const Params& a, const OrderingOptions& opt, Rng& rng)
{
    const Vertex n = g.num_vertices();
    const auto un = static_cast<std::size_t>(n);
    Scores s;
    if (opt.score == ScoreRule::Low || opt.score == ScoreRule::High)
        s = compute_scores(g, a, opt);
    std::vector<double> rnd_out;
    std::vector<double> rnd_in;
    if (opt.score == ScoreRule::Random)
    {
        rnd_out.resize(un);
        rnd_in.resize(un);
        for (auto& x : rnd_out)
            x = rng.uniform01();
        for (auto& x : rnd_in)
            x = rng.uniform01();
    }
    auto score_key = [&](Vertex v, bool out) -> double {
        const auto uv = static_cast<std::size_t>(v);
        switch (opt.score)
        {
        case ScoreRule::Low:
            return out ? s.out[uv] : s.in[uv];
        case ScoreRule::High:
            return out ? -s.out[uv] : -s.in[uv];
        case ScoreRule::Random:
            return out ? rnd_out[uv] : rnd_in[uv];
        case ScoreRule::None:
            return 0.0;
        }
        return 0.0;
    };

    std::vector<Key> kout(un);
    std::vector<Key> kin(un);
    for (Vertex v = 0; v < n; ++v)
    {
        const auto od = g.outdegree(v);
        const auto id = g.indegree(v);
        // out-order: rank desc, sinks last, indegree-1 first, then score
        {
            long long k1 = opt.use_rank ? -static_cast<long long>(scc.rank_out_of(v)) : 0;
            int k2 = 0;
            int k3 = 0;
            double k4 = score_key(v, true);
            if (opt.use_degree_rules)
            {
                if (od == 0)
                {
                    k2 = 1;
                    k4 = 0;
                }
                else if (id == 1)
                    k4 = 0;
                else
                    k3 = 1;
            }
            kout[static_cast<std::size_t>(v)] = {k1, k2, k3, k4, v};
        }
        // in-order: 2018 uses ascending in-rank (symmetric_in=false)
        {
            long long r = scc.rank_in_of(v);
            long long k1 = opt.use_rank ? (opt.symmetric_in ? -r : r) : 0;
            int k2 = 0;
            int k3 = 0;
            double k4 = score_key(v, false);
            if (opt.use_degree_rules)
            {
                if (id == 0)
                {
                    k2 = 1;
                    k4 = 0;
                }
                else if (od == 1)
                    k4 = 0;
                else
                    k3 = 1;
            }
            kin[static_cast<std::size_t>(v)] = {k1, k2, k3, k4, v};
        }
    }
    SearchOrders so;
    so.out_order.resize(un);
    so.in_order.resize(un);
    std::iota(so.out_order.begin(), so.out_order.end(), 0);
    std::iota(so.in_order.begin(), so.in_order.end(), 0);
    std::ranges::sort(so.out_order, [&](Vertex x, Vertex y) { return kout[static_cast<std::size_t>(x)] < kout[static_cast<std::size_t>(y)]; });
    std::ranges::sort(so.in_order, [&](Vertex x, Vertex y) { return kin[static_cast<std::size_t>(x)] < kin[static_cast<std::size_t>(y)]; });
    std::vector<std::int32_t> pos_out(un);
    std::vector<std::int32_t> pos_in(un);
    for (std::size_t i = 0; i < un; ++i)
    {
        pos_out[static_cast<std::size_t>(so.out_order[i])] = static_cast<std::int32_t>(i);
        pos_in[static_cast<std::size_t>(so.in_order[i])] = static_cast<std::int32_t>(i);
    }
    so.out_adj = sorted_adjacency(n, pos_out, true, g, opt.heavy_arcs_first);
    so.in_adj = sorted_adjacency(n, pos_in, false, g, opt.heavy_arcs_first);
    return so;
}

OrderedAdjacency plain_out_adjacency(const Digraph& g) { return sorted_adjacency(g.num_vertices(), {}, true, g); }
OrderedAdjacency plain_in_adjacency(const Digraph& g) { return sorted_adjacency(g.num_vertices(), {}, false, g); }

Params params_2018(int i, Rng& rng)
{
    if (i == 0)
        return {1, 4, 16, 64, 1, 4, 16, 64};
    if (i == 1)
        return {-43, 31, 11, 58, -4, 23, 43, 45};
    Params p{};
    for (auto& x : p)
        x = static_cast<double>(rng.between(-10, 89));
    return p;
}

} // namespace lsp
