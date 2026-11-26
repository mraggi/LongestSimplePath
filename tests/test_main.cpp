#include "lsp/io.hpp"
#include "lsp/solver.hpp"

#include <gtest/gtest.h>

#include <functional>
#include <set>
#include <sstream>

using namespace lsp;

namespace {

Digraph random_digraph(Vertex n, double p, Rng& rng, Weight maxw = 1)
{
    std::vector<Edge> es;
    for (Vertex u = 0; u < n; ++u)
        for (Vertex v = 0; v < n; ++v)
            if (u != v && rng.uniform01() < p)
                es.push_back({u, v, maxw == 1 ? 1 : rng.between(1, maxw)});
    return {n, es};
}

// exhaustive longest simple path (optionally from a fixed source)
Weight brute_force(const Digraph& g, Vertex source = kNoVertex)
{
    const Vertex n = g.num_vertices();
    std::vector<char> used(static_cast<std::size_t>(n), 0);
    Weight best = 0;
    std::function<void(Vertex, Weight)> rec = [&](Vertex u, Weight val) {
        best = std::max(best, val);
        for (const auto& a : g.out(u))
            if (!used[static_cast<std::size_t>(a.v)])
            {
                used[static_cast<std::size_t>(a.v)] = 1;
                rec(a.v, val + a.w);
                used[static_cast<std::size_t>(a.v)] = 0;
            }
    };
    for (Vertex s = 0; s < n; ++s)
    {
        if (source != kNoVertex && s != source)
            continue;
        used[static_cast<std::size_t>(s)] = 1;
        rec(s, 0);
        used[static_cast<std::size_t>(s)] = 0;
    }
    return best;
}

std::vector<std::vector<char>> reachability(const Digraph& g)
{
    const auto n = static_cast<std::size_t>(g.num_vertices());
    std::vector<std::vector<char>> R(n, std::vector<char>(n, 0));
    for (std::size_t s = 0; s < n; ++s)
        R[s] = reachable_from(g, static_cast<Vertex>(s));
    return R;
}

} // namespace

TEST(Graph, SelfLoopsAndParallelArcs)
{
    Digraph g(3, {{0, 1, 5}, {0, 1, 7}, {1, 1, 3}, {1, 2, 1}});
    EXPECT_EQ(g.num_edges(), 2u);
    EXPECT_EQ(*g.arc_weight(0, 1), 7);
    EXPECT_FALSE(g.arc_weight(1, 1).has_value());
    EXPECT_EQ(g.indegree(1), 1u);
    EXPECT_EQ(g.in(2)[0].v, 1);
    EXPECT_FALSE(g.is_unit_weight());
}

TEST(Graph, CheckerRejectsBadPaths)
{
    Digraph g(4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 3}, {3, 0, 4}});
    EXPECT_TRUE(check_path(g, {0, 1, 2, 3}, 6).ok);
    EXPECT_FALSE(check_path(g, {0, 1, 2, 3}, 7).ok);     // wrong value
    EXPECT_FALSE(check_path(g, {0, 2}).ok);               // missing arc
    EXPECT_FALSE(check_path(g, {0, 1, 2, 3, 0}).ok);      // repeated vertex
    EXPECT_FALSE(check_path(g, {1, 2}, std::nullopt, 0).ok); // wrong source
}

TEST(IO, EdgeListRoundTripAndAdjList)
{
    Digraph g(4, {{0, 1, 2}, {1, 2, 3}, {2, 0, 4}, {2, 3, 1}});
    std::stringstream ss;
    write_edgelist(ss, g);
    auto h = read_graph(ss, GraphFormat::EdgeList);
    EXPECT_EQ(h.num_vertices(), 4);
    EXPECT_EQ(h.num_edges(), 4u);
    EXPECT_EQ(*h.arc_weight(2, 0), 4);
    std::stringstream adj("2 3 4 5\n1 5\n1 4\n1 3\n1 2\n"); // GoGameGraphs graph1
    auto a = read_graph(adj, GraphFormat::AdjList);
    EXPECT_EQ(a.num_vertices(), 5);
    EXPECT_EQ(a.num_edges(), 12u);
    EXPECT_TRUE(a.arc_weight(0, 4).has_value());
}

TEST(Structure, SccMatchesReachabilityAndRanksAreLongestCondensationPaths)
{
    Rng rng(42);
    for (int t = 0; t < 200; ++t)
    {
        const auto n = static_cast<Vertex>(rng.between(1, 25));
        auto g = random_digraph(n, rng.uniform01() * 0.2, rng);
        auto scc = compute_scc(g);
        auto R = reachability(g);
        for (Vertex u = 0; u < n; ++u)
            for (Vertex v = 0; v < n; ++v)
            {
                const bool same = R[u][v] && R[v][u];
                EXPECT_EQ(same, scc.comp[u] == scc.comp[v]);
            }
        // topological ids
        for (Vertex u = 0; u < n; ++u)
            for (const auto& a : g.out(u))
                EXPECT_LE(scc.comp[u], scc.comp[a.v]);
        // blocks
        for (std::int32_t c = 0; c < scc.num_components; ++c)
            for (auto i = scc.start[c]; i < scc.start[c + 1]; ++i)
                EXPECT_EQ(scc.comp[scc.block[i]], c);
        // ranks by brute force over the condensation DAG
        const auto C = scc.num_components;
        std::vector<std::set<std::int32_t>> cadj(C);
        for (Vertex u = 0; u < n; ++u)
            for (const auto& a : g.out(u))
                if (scc.comp[u] != scc.comp[a.v])
                    cadj[scc.comp[u]].insert(scc.comp[a.v]);
        std::function<int(int)> lo = [&](int c) {
            int b = 0;
            for (int d : cadj[c])
                b = std::max(b, 1 + lo(d));
            return b;
        };
        for (int c = 0; c < C; ++c)
            EXPECT_EQ(scc.rank_out[c], lo(c));
    }
}

TEST(PTO, DpEqualsLongestPathOfForwardSubgraph)
{
    Rng rng(7);
    for (int t = 0; t < 300; ++t)
    {
        const auto n = static_cast<Vertex>(rng.between(1, 9));
        auto g = random_digraph(n, 0.35, rng, 9);
        auto scc = compute_scc(g);
        auto so = make_orders(g, scc, params_2018(0, rng), {}, rng);
        PseudoTopoOrder T(g, scc, so.in_adj, so.out_adj);
        T.randomize(rng);
        // forward subgraph
        std::vector<Edge> fe;
        for (Vertex u = 0; u < n; ++u)
            for (const auto& a : g.out(u))
                if (T.position(u) < T.position(a.v))
                    fe.push_back({u, a.v, a.w});
        Digraph f(n, fe);
        const Weight bf = brute_force(f);
        EXPECT_EQ(T.value(), bf);
        auto p = T.best_path();
        auto chk = check_path(f, p.vertices, p.value);
        EXPECT_TRUE(chk.ok) << chk.message;
        EXPECT_TRUE(T.is_strong_pto());
    }
}

TEST(PTO, ImpositionKeepsStrongPtoAndPathIsForward)
{
    Rng rng(11);
    for (int t = 0; t < 300; ++t)
    {
        const auto n = static_cast<Vertex>(rng.between(2, 40));
        auto g = random_digraph(n, 3.0 / n, rng, 5);
        auto scc = compute_scc(g);
        auto so = make_orders(g, scc, params_2018(0, rng), {}, rng);
        // a random simple path via a randomized greedy walk
        std::vector<Vertex> p{static_cast<Vertex>(rng.below(n))};
        std::vector<char> used(n, 0);
        used[p[0]] = 1;
        while (true)
        {
            std::vector<Vertex> c;
            for (const auto& a : g.out(p.back()))
                if (!used[a.v])
                    c.push_back(a.v);
            if (c.empty())
                break;
            p.push_back(c[rng.below(c.size())]);
            used[p.back()] = 1;
        }
        const Weight w = check_path(g, p).recomputed;
        PseudoTopoOrder T(g, scc, so.in_adj, so.out_adj);
        T.randomize(rng);
        T.impose(p);
        EXPECT_TRUE(T.is_strong_pto());
        for (std::size_t i = 1; i < p.size(); ++i)
            EXPECT_LT(T.position(p[i - 1]), T.position(p[i]));
        EXPECT_GE(T.value(), w);
        // opening passes never decrease the value and keep the order a strong PTO
        OpenerOptions oo;
        oo.skip_block_size = 0;
        oo.min_free = 1;
        oo.sort_iters = 50;
        OpenStats st;
        if (t % 2)
            oo.border = BOrderRule::DfsPostorder;
        for (int k = 0; k < 5; ++k)
        {
            const Weight before = T.value();
            T.open_pass(oo, rng, Deadline::never(), st);
            EXPECT_GE(T.value(), before);
            EXPECT_TRUE(T.is_strong_pto());
            auto bp = T.best_path();
            EXPECT_TRUE(check_path(g, bp.vertices, bp.value).ok);
        }
    }
}

TEST(PTO, DfsReversePostorderIsWeakPtoOfInducedSubgraph)
{
    Rng rng(17);
    for (int t = 0; t < 200; ++t)
    {
        const auto n = static_cast<Vertex>(rng.between(2, 60));
        auto g = random_digraph(n, 2.5 / n, rng);
        auto scc = compute_scc(g);
        auto so = make_orders(g, scc, params_2018(0, rng), {}, rng);
        PseudoTopoOrder T(g, scc, so.in_adj, so.out_adj);
        T.randomize(rng);
        // a random window of positions (a set B that is generally not strongly connected)
        const auto b = static_cast<std::int32_t>(rng.below(n));
        const auto c = b + 1 + static_cast<std::int32_t>(rng.below(n - b));
        T.order_free_dfs(b, c, b > 0 ? T.at(b - 1) : kNoVertex, rng);
        std::vector<char> keep(n, 0);
        for (auto i = b; i < c; ++i)
            keep[T.at(i)] = 1;
        std::vector<Vertex> old;
        auto h = g.induced(keep, old);
        auto hs = compute_scc(h);
        std::vector<std::int32_t> newid(n, -1);
        for (std::size_t k = 0; k < old.size(); ++k)
            newid[old[k]] = static_cast<std::int32_t>(k);
        for (Vertex u = 0; u < n; ++u)
            for (const auto& a : g.out(u))
                if (keep[u] && keep[a.v] && T.position(a.v) < T.position(u))
                    EXPECT_EQ(hs.comp[newid[u]], hs.comp[newid[a.v]]) << "backward arc between SCCs of D[B]";
        for (auto i = b; i < c; ++i)
            EXPECT_EQ(T.position(T.at(i)), i);
        // Prop. (c): every vertex of B reachable from the anchor through B is the end of a
        // forward path anchor -> r1 -> ... -> r inside B
        if (b > 0)
        {
            const Vertex anchor = T.at(b - 1);
            std::vector<char> reach_fwd(n, 0), reach_any(n, 0);
            std::vector<Vertex> st{anchor};
            reach_any[anchor] = 1;
            while (!st.empty())
            {
                const Vertex u = st.back();
                st.pop_back();
                for (const auto& a : g.out(u))
                    if (keep[a.v] && !reach_any[a.v])
                    {
                        reach_any[a.v] = 1;
                        st.push_back(a.v);
                    }
            }
            // forward reachability: process positions in increasing order
            for (const auto& a : g.out(anchor))
                if (keep[a.v])
                    reach_fwd[a.v] = 1;
            for (auto i = b; i < c; ++i)
            {
                const Vertex u = T.at(i);
                if (!reach_fwd[u])
                    continue;
                for (const auto& a : g.out(u))
                    if (keep[a.v] && T.position(a.v) > i)
                        reach_fwd[a.v] = 1;
            }
            for (auto i = b; i < c; ++i)
                if (reach_any[T.at(i)])
                    EXPECT_TRUE(reach_fwd[T.at(i)]) << "reachable free vertex without a forward path from the anchor";
        }
    }
}

TEST(PTO, IncrementalDpEqualsFullDp)
{
    Rng rng(3);
    for (int t = 0; t < 50; ++t)
    {
        const auto n = static_cast<Vertex>(rng.between(5, 60));
        auto g = random_digraph(n, 4.0 / n, rng, 7);
        auto scc = compute_scc(g);
        auto so = make_orders(g, scc, params_2018(0, rng), {}, rng);
        PseudoTopoOrder A(g, scc, so.in_adj, so.out_adj);
        PseudoTopoOrder B(g, scc, so.in_adj, so.out_adj);
        B.incremental = false;
        A.randomize(rng);
        for (int k = 0; k < 100; ++k)
        {
            // random swap inside a block
            const auto c = static_cast<std::int32_t>(rng.below(scc.num_components));
            if (scc.size(c) < 2)
                continue;
            const auto i = scc.start[c] + static_cast<std::int32_t>(rng.below(scc.size(c)));
            const auto j = scc.start[c] + static_cast<std::int32_t>(rng.below(scc.size(c)));
            A.swap_positions(i, j);
            B.impose(A.order()); // impose the whole order
            EXPECT_EQ(A.value(), B.value());
        }
    }
}

TEST(DFS, ExhaustiveEnumerationFindsBestRootedPath)
{
    Rng rng(5);
    for (int t = 0; t < 200; ++t)
    {
        const auto n = static_cast<Vertex>(rng.between(1, 9));
        auto g = random_digraph(n, 0.4, rng, 6);
        auto adj = plain_out_adjacency(g);
        const auto s = static_cast<Vertex>(rng.below(n));
        InPlaceDfs d(n, adj);
        d.reset({s});
        d.greedy_extend();
        DfsStats st;
        auto best = d.run({}, Deadline::never(), st, nullptr);
        auto chk = check_path(g, best);
        EXPECT_TRUE(chk.ok) << chk.message;
        EXPECT_EQ(best.front(), s);
        EXPECT_EQ(chk.recomputed, brute_force(g, s));
    }
}

TEST(Solver, AllMethodsReturnValidPathsAndClassicIsOptimalOnTinyGraphs)
{
    Rng rng(9);
    int optimal = 0;
    for (int t = 0; t < 40; ++t)
    {
        const auto n = static_cast<Vertex>(rng.between(2, 9));
        auto g = random_digraph(n, 0.35, rng, t % 2 ? 1 : 9);
        const Weight opt = brute_force(g);
        for (const std::string m : {"classic", "ils", "rdfs", "greedy", "gamm"})
        {
            SolverOptions o;
            o.apply_preset("2018");
            o.method = m;
            o.seed = static_cast<std::uint64_t>(t);
            o.largest_wcc_only = false;
            o.dfs_stop = {0, 2000};
            o.pto_stop = {0, 3};
            o.param_restarts = 2;
            o.time_limit = (m == "classic") ? 0 : 0.05;
            auto r = solve(g, o);
            auto chk = check_path(g, r.best.vertices, r.best.value);
            EXPECT_TRUE(chk.ok) << m << ": " << chk.message;
            EXPECT_LE(r.best.value, opt);
            ASSERT_FALSE(r.trajectory.empty());
            EXPECT_EQ(r.trajectory.back().value, r.best.value) << m;
            if (m == "classic" && r.best.value == opt)
                ++optimal;
        }
    }
    // Not guaranteed (the DFS continues an enumeration instead of restarting it), but it
    // should almost always be the case on graphs with at most 9 vertices.
    EXPECT_GE(optimal, 36);
}

TEST(Solver, RootedVariantStartsAtSource)
{
    Rng rng(21);
    for (int t = 0; t < 40; ++t)
    {
        const auto n = static_cast<Vertex>(rng.between(2, 9));
        auto g = random_digraph(n, 0.4, rng);
        const auto s = static_cast<Vertex>(rng.below(n));
        SolverOptions o;
        o.apply_preset("2018");
        o.time_limit = 0;
        o.source = s;
        o.dfs_stop = {0, 5000};
        o.pto_stop = {0, 3};
        o.param_restarts = 2;
        auto r = solve(g, o);
        auto chk = check_path(g, r.best.vertices, r.best.value, s);
        EXPECT_TRUE(chk.ok) << chk.message;
        EXPECT_EQ(r.best.value, brute_force(g, s));
    }
}

TEST(Solver, WorkBasedRunsAreDeterministic)
{
    Rng rng(1);
    auto g = random_digraph(300, 0.01, rng, 20);
    SolverOptions o;
    o.apply_preset("2018");
    o.time_limit = 0;
    o.seed = 77;
    o.dfs_stop = {0, 3000};
    o.pto_stop = {0, 2};
    o.param_restarts = 3;
    auto a = solve(g, o);
    auto b = solve(g, o);
    EXPECT_EQ(a.best.vertices, b.best.vertices);
    EXPECT_EQ(a.best.value, b.best.value);
    EXPECT_TRUE(check_path(g, a.best.vertices, a.best.value).ok);
}

TEST(Options, PresetsKeepTimeSeedAndSource)
{
    for (const std::string p : {"2018", "new", "newap"})
    {
        SolverOptions o;
        o.set("time", "7.5");
        o.set("seed", "42");
        o.set("source", "3");
        o.set("border", "shuffle");
        o.set("preset", p);
        EXPECT_EQ(o.time_limit, 7.5) << p;
        EXPECT_EQ(o.seed, 42u) << p;
        EXPECT_EQ(o.source, 3) << p;
        EXPECT_EQ(o.opener.border, p != "2018" ? BOrderRule::DfsPostorder : BOrderRule::HeuristicFirst) << p;
    }
}

TEST(Options, DefaultIsNewapWithTimeLimit)
{
    const SolverOptions d;
    SolverOptions n;
    n.apply_preset("newap");
    n.time_limit = SolverOptions::kDefaultTimeLimit;
    EXPECT_EQ(d.describe(), n.describe());
    EXPECT_EQ(d.method, "ils");
    EXPECT_TRUE(d.ap_init);
    EXPECT_FALSE(d.largest_wcc_only);
    EXPECT_GT(d.time_limit, 0);
    SolverOptions o;
    o.apply_preset("2018"); // the 2016-2018 method is still available
    EXPECT_EQ(o.method, "classic");
    EXPECT_FALSE(o.ap_init);
    EXPECT_TRUE(o.largest_wcc_only);
    EXPECT_EQ(o.opener.border, BOrderRule::HeuristicFirst);
    EXPECT_EQ(o.time_limit, SolverOptions::kDefaultTimeLimit); // presets keep the time limit
}

TEST(Ordering, HeavyFirstDoesNotChangeUnitWeightRuns)
{
    Rng rng(5);
    auto g = random_digraph(300, 0.012, rng);
    SolverOptions a;
    a.apply_preset("new");
    a.dfs_stop = {0, 3000};
    a.pto_stop = {0, 2};
    a.time_limit = 0;
    a.method = "classic";
    a.param_restarts = 2;
    SolverOptions b = a;
    b.ordering.heavy_arcs_first = false;
    EXPECT_EQ(solve(g, a).best.vertices, solve(g, b).best.vertices);
}

namespace {

// exhaustive AP optimum: permutations of V+D (D = n) with the allowed choices of ap.hpp
Weight brute_force_ap(const Digraph& g, Vertex source)
{
    const Vertex n = g.num_vertices();
    const bool rooted = source != kNoVertex;
    std::vector<char> used(static_cast<std::size_t>(n) + 1, 0);
    Weight best = kNegInf;
    std::function<void(Vertex, Weight)> rec = [&](Vertex u, Weight val) {
        if (u > n)
        {
            best = std::max(best, val);
            return;
        }
        for (Vertex v = 0; v <= n; ++v)
        {
            if (used[static_cast<std::size_t>(v)])
                continue;
            Weight w = 0;
            if (u == n) // D chooses the first vertex (or D: empty path)
            {
                if (rooted ? v != source : false)
                    continue;
            }
            else if (v == n) // u ends the path
            {
            }
            else if (v == u) // u skipped
            {
                if (rooted && u == source)
                    continue;
            }
            else
            {
                auto aw = g.arc_weight(u, v);
                if (!aw || (rooted && v == source))
                    continue;
                w = *aw;
            }
            used[static_cast<std::size_t>(v)] = 1;
            rec(u + 1, val + w);
            used[static_cast<std::size_t>(v)] = 0;
        }
    };
    rec(0, 0);
    return best;
}

Digraph random_signed_digraph(Vertex n, double p, Rng& rng, Weight lo, Weight hi)
{
    std::vector<Edge> es;
    for (Vertex u = 0; u < n; ++u)
        for (Vertex v = 0; v < n; ++v)
            if (u != v && rng.uniform01() < p)
                es.push_back({u, v, rng.between(lo, hi)});
    return {n, es};
}

} // namespace

TEST(AP, BoundEqualsExhaustiveApOptimumAndBoundsLongestPath)
{
    Rng rng(31);
    int strict = 0;
    int negative = 0;
    for (int t = 0; t < 300; ++t)
    {
        const auto n = static_cast<Vertex>(rng.between(1, 7));
        const int kind = t % 3; // unit, positive, signed weights
        auto g = kind == 0 ? random_digraph(n, 0.4, rng) : kind == 1 ? random_digraph(n, 0.4, rng, 9)
                                                                      : random_signed_digraph(n, 0.45, rng, -6, 9);
        const Vertex src = (t % 2) ? static_cast<Vertex>(rng.below(static_cast<std::uint64_t>(n))) : kNoVertex;
        auto ap = solve_ap(g, src);
        ASSERT_TRUE(ap.solved);
        EXPECT_TRUE(ap.certified) << t;
        EXPECT_EQ(ap.bound, brute_force_ap(g, src)) << t;
        const Weight lp = brute_force(g, src);
        EXPECT_GE(ap.bound, lp) << t;
        strict += ap.bound > lp ? 1 : 0;
        for (const auto& e : g.edges())
            negative += e.weight < 0 ? 1 : 0;
        // the AP solution is consistent: path + cycles are disjoint, use arcs, and add up to the bound
        std::vector<char> seen(static_cast<std::size_t>(n), 0);
        Weight total = 0;
        auto walk = [&](const std::vector<Vertex>& vs, bool closed) {
            for (std::size_t i = 0; i < vs.size(); ++i)
            {
                EXPECT_FALSE(seen[static_cast<std::size_t>(vs[i])]);
                seen[static_cast<std::size_t>(vs[i])] = 1;
                if (i + 1 < vs.size() || (closed && vs.size() > 1))
                {
                    auto w = g.arc_weight(vs[i], vs[(i + 1) % vs.size()]);
                    ASSERT_TRUE(w.has_value());
                    total += *w;
                }
            }
        };
        walk(ap.path, false);
        for (const auto& c : ap.cycles)
            walk(c, true);
        EXPECT_EQ(total, ap.bound) << t;
        if (src != kNoVertex)
        {
            ASSERT_FALSE(ap.path.empty());
            EXPECT_EQ(ap.path.front(), src);
        }
        PatchStats ps;
        auto p = patch_ap(g, ap, src, &ps);
        if (p.empty())
            continue;
        auto chk = check_path(g, p, std::nullopt, src != kNoVertex ? std::optional<Vertex>(src) : std::nullopt);
        EXPECT_TRUE(chk.ok) << t << ": " << chk.message;
        EXPECT_LE(chk.recomputed, lp);
    }
    // the test must exercise relaxations that are not tight, and negative arcs
    EXPECT_GT(strict, 30);
    EXPECT_GT(negative, 100);
    std::cout << "[ info ] AP bound > longest path in " << strict << "/300 cases\n";
}

TEST(AP, PatchingMergesCyclesIntoOnePathOnLargerGraphs)
{
    Rng rng(8);
    for (int t = 0; t < 20; ++t)
    {
        auto g = random_digraph(400, 0.008, rng, t % 2 ? 1 : 50);
        const Vertex src = t % 4 == 3 ? 0 : kNoVertex;
        auto ap = solve_ap(g, src);
        ASSERT_TRUE(ap.solved);
        EXPECT_TRUE(ap.certified);
        PatchStats ps;
        auto p = patch_ap(g, ap, src, &ps);
        auto chk = check_path(g, p, std::nullopt, src != kNoVertex ? std::optional<Vertex>(src) : std::nullopt);
        EXPECT_TRUE(chk.ok) << chk.message;
        EXPECT_LE(chk.recomputed, ap.bound);
        EXPECT_EQ(static_cast<std::size_t>(ps.merged + ps.lossy + ps.left), ap.cycles.size());
    }
}

TEST(Solver, NewApPresetIsValidAndProvenOptimalOnlyWhenOptimal)
{
    Rng rng(17);
    int proven = 0;
    for (int t = 0; t < 60; ++t)
    {
        const auto n = static_cast<Vertex>(rng.between(2, 9));
        auto g = t % 3 == 2 ? random_signed_digraph(n, 0.4, rng, -5, 9) : random_digraph(n, 0.35, rng, t % 2 ? 1 : 9);
        const Vertex src = t % 4 == 1 ? static_cast<Vertex>(rng.below(static_cast<std::uint64_t>(n))) : kNoVertex;
        SolverOptions o;
        o.set("preset", "newap");
        o.seed = static_cast<std::uint64_t>(t);
        o.source = src;
        o.time_limit = 0.05;
        auto r = solve(g, o);
        auto chk = check_path(g, r.best.vertices, r.best.value, src != kNoVertex ? std::optional<Vertex>(src) : std::nullopt);
        EXPECT_TRUE(chk.ok) << chk.message;
        const Weight opt = brute_force(g, src);
        EXPECT_LE(r.best.value, opt);
        EXPECT_TRUE(r.ap.ran && r.ap.solved && r.ap.certified && r.ap.bound_global);
        EXPECT_GE(r.ap.bound, opt);
        if (r.proven_optimal)
        {
            EXPECT_EQ(r.best.value, opt);
            ++proven;
        }
    }
    EXPECT_GT(proven, 10);
}
