#include "lsp/solver.hpp"

#include <fstream>

#include <sstream>
#include <stdexcept>

namespace lsp {

// ---------------------------------------------------------------- options

namespace {

bool parse_bool(const std::string& v)
{
    if (v == "1" || v == "true" || v == "on" || v == "yes")
        return true;
    if (v == "0" || v == "false" || v == "off" || v == "no")
        return false;
    throw std::invalid_argument("bad boolean: " + v);
}

ScoreRule parse_score(const std::string& v)
{
    if (v == "low")
        return ScoreRule::Low;
    if (v == "high")
        return ScoreRule::High;
    if (v == "random")
        return ScoreRule::Random;
    if (v == "none")
        return ScoreRule::None;
    throw std::invalid_argument("bad score rule: " + v);
}

std::string score_name(ScoreRule r)
{
    switch (r)
    {
    case ScoreRule::Low:
        return "low";
    case ScoreRule::High:
        return "high";
    case ScoreRule::Random:
        return "random";
    case ScoreRule::None:
        return "none";
    }
    return "?";
}

BOrderRule parse_border(const std::string& v)
{
    if (v == "first")
        return BOrderRule::HeuristicFirst;
    if (v == "uniform")
        return BOrderRule::HeuristicUniform;
    if (v == "shuffle")
        return BOrderRule::ShuffleOnly;
    if (v == "dfs")
        return BOrderRule::DfsPostorder;
    throw std::invalid_argument("bad B-order rule: " + v);
}

std::string border_name(BOrderRule r)
{
    switch (r)
    {
    case BOrderRule::HeuristicFirst:
        return "first";
    case BOrderRule::HeuristicUniform:
        return "uniform";
    case BOrderRule::ShuffleOnly:
        return "shuffle";
    case BOrderRule::DfsPostorder:
        return "dfs";
    }
    return "?";
}

template <class T>
std::string str(const T& x)
{
    std::ostringstream ss;
    ss << x;
    return ss.str();
}

} // namespace

void SolverOptions::apply_preset(const std::string& name)
{
    if (name != "2018" && name != "new" && name != "newap")
        throw std::invalid_argument("unknown preset: " + name);
    // a preset resets the algorithmic options but keeps the run parameters
    const double t = time_limit;
    const std::uint64_t s = seed;
    const Vertex src = source;
    *this = SolverOptions{};
    time_limit = t;
    seed = s;
    source = src;
    if (name == "2018")
        return;
    method = "ils";
    largest_wcc_only = false;
    ordering.symmetric_in = true;
    ordering.in_score_typo = false;
    ordering.heavy_arcs_first = true; // weighted graphs: try heavier arcs first (no effect with unit weights)
    faithful_backward_greedy = false;
    dfs_stop = {0, 100000};
    pto_stop = {0, 3};
    opener.border = BOrderRule::DfsPostorder;
    opener.skip_block_size = 1; // open every non-trivial SCC (needed for the certificates)
    opener.min_free = 1;
    ap_init = name == "newap";
}

void SolverOptions::set(const std::string& key, const std::string& value)
{
    auto d = [&] { return std::stod(value); };
    auto i64 = [&] { return std::stoll(value); };
    auto b = [&] { return parse_bool(value); };
    if (key == "preset")
        apply_preset(value);
    else if (key == "method")
        method = value;
    else if (key == "time")
        time_limit = d();
    else if (key == "seed")
        seed = std::stoull(value);
    else if (key == "source")
        source = static_cast<Vertex>(i64());
    else if (key == "largest_wcc")
        largest_wcc_only = b();
    else if (key == "rank")
        ordering.use_rank = b();
    else if (key == "degree_rules")
        ordering.use_degree_rules = b();
    else if (key == "score")
        ordering.score = parse_score(value);
    else if (key == "symmetric_in")
        ordering.symmetric_in = b();
    else if (key == "in_score_typo")
        ordering.in_score_typo = b();
    else if (key == "weighted_scores")
        ordering.weighted_scores = b();
    else if (key == "heavy_first")
        ordering.heavy_arcs_first = b();
    else if (key == "param_restarts")
        param_restarts = static_cast<int>(i64());
    else if (key == "params")
        params = value;
    else if (key == "fwd_starts")
        fwd_starts = static_cast<int>(i64());
    else if (key == "bwd_starts")
        bwd_starts = static_cast<int>(i64());
    else if (key == "erase")
        erase = static_cast<int>(i64());
    else if (key == "backward")
        backward = b();
    else if (key == "other_direction")
        other_direction = b();
    else if (key == "ranked_starts")
        ranked_starts = b();
    else if (key == "faithful_backward_greedy")
        faithful_backward_greedy = b();
    else if (key == "dfs_seconds")
        dfs_stop.seconds_without_improvement = d();
    else if (key == "dfs_steps")
        dfs_stop.steps_without_improvement = i64();
    else if (key == "pto")
        pto = b();
    else if (key == "pto_each_dfs")
        pto_each_dfs = b();
    else if (key == "pto_mode")
    {
        if (value == "open")
            pto_mode = PtoMode::Open;
        else if (value == "reimpose")
            pto_mode = PtoMode::Reimpose;
        else
            throw std::invalid_argument("bad pto_mode: " + value);
    }
    else if (key == "pto_seconds")
        pto_stop.seconds_without_improvement = d();
    else if (key == "pto_passes")
        pto_stop.steps_without_improvement = i64();
    else if (key == "pto_incremental")
        pto_incremental = b();

    else if (key == "skip_block")
        opener.skip_block_size = static_cast<std::int32_t>(i64());
    else if (key == "min_free")
        opener.min_free = static_cast<std::int32_t>(i64());
    else if (key == "border")
        opener.border = parse_border(value);
    else if (key == "sort_iters")
        opener.sort_iters = i64();
    else if (key == "sort_iters_per_free")
        opener.sort_iters_per_free = d();
    else if (key == "sweep")
        opener.sweep = b();
    else if (key == "random_gaps")
        opener.random_gap_order = b();
    else if (key == "reshuffle_each_gap")
        opener.reshuffle_each_gap = b();
    else if (key == "try_reverse")
        opener.try_reverse = b();
    else if (key == "free_first")
        opener.free_first = b();
    else if (key == "sideways")
        opener.sideways = d();
    else if (key == "sideways_max")
        sideways_max = i64();
    else if (key == "ils_keep")
        ils_keep_fraction = d();
    else if (key == "ils_dfs_steps")
        ils_dfs_steps = i64();
    else if (key == "ils_rounds")
        ils_rounds = i64();
    else if (key == "init_path")
        init_path = value;
    else if (key == "ap_init")
        ap_init = b();
    else if (key == "rdfs_steps")
        rdfs_steps = i64();
    else if (key == "ga_population")
        ga_population = static_cast<int>(i64());
    else if (key == "ga_elite")
        ga_elite = d();
    else if (key == "ga_stagnation")
        ga_stagnation = static_cast<int>(i64());
    else
        throw std::invalid_argument("unknown option: " + key);
}

std::map<std::string, std::string> SolverOptions::describe() const
{
    std::map<std::string, std::string> m;
    m["method"] = method;
    m["time"] = str(time_limit);
    m["seed"] = str(seed);
    m["source"] = str(source);
    m["largest_wcc"] = str(largest_wcc_only);
    m["rank"] = str(ordering.use_rank);
    m["degree_rules"] = str(ordering.use_degree_rules);
    m["score"] = score_name(ordering.score);
    m["symmetric_in"] = str(ordering.symmetric_in);
    m["in_score_typo"] = str(ordering.in_score_typo);
    m["weighted_scores"] = str(ordering.weighted_scores);
    m["heavy_first"] = str(ordering.heavy_arcs_first);
    m["param_restarts"] = str(param_restarts);
    m["params"] = params;
    m["fwd_starts"] = str(fwd_starts);
    m["bwd_starts"] = str(bwd_starts);
    m["erase"] = str(erase);
    m["backward"] = str(backward);
    m["other_direction"] = str(other_direction);
    m["ranked_starts"] = str(ranked_starts);
    m["faithful_backward_greedy"] = str(faithful_backward_greedy);
    m["dfs_seconds"] = str(dfs_stop.seconds_without_improvement);
    m["dfs_steps"] = str(dfs_stop.steps_without_improvement);
    m["pto"] = str(pto);
    m["pto_each_dfs"] = str(pto_each_dfs);
    m["pto_mode"] = pto_mode == PtoMode::Open ? "open" : "reimpose";
    m["pto_seconds"] = str(pto_stop.seconds_without_improvement);
    m["pto_passes"] = str(pto_stop.steps_without_improvement);
    m["pto_incremental"] = str(pto_incremental);
    m["skip_block"] = str(opener.skip_block_size);
    m["min_free"] = str(opener.min_free);
    m["border"] = border_name(opener.border);
    m["sort_iters"] = str(opener.sort_iters);
    m["sort_iters_per_free"] = str(opener.sort_iters_per_free);
    m["sweep"] = str(opener.sweep);
    m["random_gaps"] = str(opener.random_gap_order);
    m["reshuffle_each_gap"] = str(opener.reshuffle_each_gap);
    m["try_reverse"] = str(opener.try_reverse);
    m["free_first"] = str(opener.free_first);
    m["sideways"] = str(opener.sideways);
    m["sideways_max"] = str(sideways_max);
    m["ils_keep"] = str(ils_keep_fraction);
    m["ils_dfs_steps"] = str(ils_dfs_steps);
    m["ils_rounds"] = str(ils_rounds);
    if (!init_path.empty())
        m["init_path"] = init_path;
    m["ap_init"] = str(ap_init);
    m["rdfs_steps"] = str(rdfs_steps);
    m["ga_population"] = str(ga_population);
    m["ga_elite"] = str(ga_elite);
    m["ga_stagnation"] = str(ga_stagnation);
    return m;
}

// ---------------------------------------------------------------- solver

namespace {

class Solver
{
public:
    Solver(const Digraph& input, const SolverOptions& opt) : opt_(opt), rng_(opt.seed)
    {
        deadline_ = Deadline::in(opt.time_limit);
        prepare(input);
        res_.preprocess_seconds = clock_.seconds();
        res_.work_vertices = g_.num_vertices();
        res_.work_edges = g_.num_edges();
    }

    SolveResult run()
    {
        if (g_.num_vertices() > 0)
        {
            // the trivial path is always available
            offer({source_ != kNoVertex ? source_ : 0}, 0);
            if (opt_.method == "classic")
                classic();
            else if (opt_.method == "ils")
                ils();
            else if (opt_.method == "rdfs")
                rdfs();
            else if (opt_.method == "greedy")
                greedy();
            else if (opt_.method == "gamm")
                gamm();
            else
                throw std::invalid_argument("unknown method: " + opt_.method);
        }
        res_.best.vertices.clear();
        for (Vertex v : best_)
            res_.best.vertices.push_back(old_of_new_[static_cast<std::size_t>(v)]);
        res_.best.value = best_value_;
        res_.proven_optimal = proven_optimal();
        res_.total_seconds = clock_.seconds();
        return std::move(res_);
    }

private:
    // ---- preparation
    void prepare(const Digraph& input)
    {
        const Vertex n = input.num_vertices();
        std::vector<char> keep(static_cast<std::size_t>(n), 1);
        if (opt_.source != kNoVertex)
        {
            if (opt_.source < 0 || opt_.source >= n)
                throw std::invalid_argument("source out of range");
            keep = reachable_from(input, opt_.source);
        }
        else if (opt_.largest_wcc_only && n > 0)
        {
            std::int32_t cnt = 0;
            auto w = weak_components(input, cnt);
            std::vector<std::int32_t> sz(static_cast<std::size_t>(cnt), 0);
            for (auto c : w)
                ++sz[static_cast<std::size_t>(c)];
            const auto big = static_cast<std::int32_t>(std::ranges::max_element(sz) - sz.begin());
            for (Vertex v = 0; v < n; ++v)
                keep[static_cast<std::size_t>(v)] = static_cast<char>(w[static_cast<std::size_t>(v)] == big);
        }
        g_ = input.induced(keep, old_of_new_);
        if (opt_.source != kNoVertex)
        {
            for (Vertex v = 0; v < g_.num_vertices(); ++v)
                if (old_of_new_[static_cast<std::size_t>(v)] == opt_.source)
                    source_ = v;
            g_ = g_.without_in_arcs(source_);
        }
        // a bound on g_ bounds the whole problem unless vertices were dropped for another reason
        // than unreachability from the source
        bound_global_ = opt_.source != kNoVertex || g_.num_vertices() == n;
        scc_ = compute_scc(g_);
    }

    // ---- incumbent
    void record(Weight v)
    {
        if (v > best_value_)
        {
            best_value_ = v;
            res_.trajectory.push_back({clock_.seconds(), v});
        }
    }
    void offer(const std::vector<Vertex>& p, Weight v)
    {
        if (p.empty())
            return;
        if (v > best_path_value_)
        {
            best_ = p;
            best_path_value_ = v;
        }
        record(v);
    }
    ValueSink sink()
    {
        return [this](Weight v) { record(v); };
    }
    Weight path_value(const std::vector<Vertex>& p) const
    {
        Weight s = 0;
        for (std::size_t i = 1; i < p.size(); ++i)
            s += *g_.arc_weight(p[i - 1], p[i]);
        return s;
    }

    // ---- DFS building blocks (2018 dfs_search_path_forward/backward)
    std::vector<Vertex> dfs_forward(const SearchOrders& so, std::vector<Vertex> start, StopRule stop)
    {
        InPlaceDfs d(g_.num_vertices(), so.out_adj);
        d.reset(start);
        d.greedy_extend();
        auto p = d.run(stop, deadline_, res_.dfs, sink());
        offer(p, path_value(p));
        return p;
    }
    std::vector<Vertex> dfs_backward(const SearchOrders& so, std::vector<Vertex> start, StopRule stop)
    {
        if (opt_.faithful_backward_greedy)
        {
            // 2018: ExpandGreedyBack before the backward enumeration
            InPlaceDfs f(g_.num_vertices(), so.out_adj);
            f.reset(start);
            f.greedy_extend();
            start = f.sequence();
        }
        InPlaceDfs d(g_.num_vertices(), so.in_adj);
        d.reset(reversed(std::move(start)));
        d.greedy_extend();
        auto p = reversed(d.run(stop, deadline_, res_.dfs, sink()));
        offer(p, path_value(p));
        return p;
    }

    Vertex random_vertex() { return static_cast<Vertex>(rng_.below(static_cast<std::uint64_t>(g_.num_vertices()))); }

    std::vector<Vertex> dfs_forward_full(const SearchOrders& so)
    {
        if (source_ != kNoVertex)
        {
            auto p = dfs_forward(so, {source_}, opt_.dfs_stop);
            offer(p, path_value(p));
            return p;
        }
        std::vector<Vertex> best;
        Weight bv = kNegInf;
        const int starts = std::min<int>(opt_.fwd_starts, g_.num_vertices());
        for (int i = 0; i < starts && !deadline_.expired(); ++i)
        {
            Vertex s = so.out_order[static_cast<std::size_t>(i / 2)];
            if (!opt_.ranked_starts || i % 2 == 1)
                s = random_vertex();
            auto p = dfs_forward(so, {s}, opt_.dfs_stop);
            const Weight v = path_value(p);
            if (v > bv)
            {
                bv = v;
                best = std::move(p);
            }
        }
        offer(best, bv);
        if (opt_.other_direction && !deadline_.expired())
        {
            auto p = best;
            int k = opt_.erase;
            if (2 * k > static_cast<int>(p.size()))
                k = static_cast<int>(p.size()) / 2;
            p.erase(p.begin(), p.begin() + k);
            p = dfs_backward(so, p, opt_.dfs_stop);
            const Weight v = path_value(p);
            if (v > bv)
            {
                bv = v;
                best = std::move(p);
                offer(best, bv);
            }
        }
        return best;
    }

    std::vector<Vertex> dfs_backward_full(const SearchOrders& so)
    {
        std::vector<Vertex> best;
        Weight bv = kNegInf;
        const int starts = std::min<int>(opt_.bwd_starts, g_.num_vertices());
        for (int i = 0; i < starts && !deadline_.expired(); ++i)
        {
            Vertex s = so.in_order[static_cast<std::size_t>(i / 2)];
            if (!opt_.ranked_starts || i % 2 == 1)
                s = random_vertex();
            auto p = dfs_backward(so, {s}, opt_.dfs_stop);
            const Weight v = path_value(p);
            if (v > bv)
            {
                bv = v;
                best = std::move(p);
            }
        }
        if (best.empty())
            return best;
        offer(best, bv);
        if (opt_.other_direction && !deadline_.expired())
        {
            auto p = best;
            int k = opt_.erase;
            if (k > static_cast<int>(p.size()) / 2)
                k = static_cast<int>(p.size()) / 2;
            p.resize(p.size() - static_cast<std::size_t>(k));
            p = dfs_forward(so, p, opt_.dfs_stop);
            const Weight v = path_value(p);
            if (v > bv)
            {
                bv = v;
                best = std::move(p);
                offer(best, bv);
            }
        }
        return best;
    }

    // ---- PTO search (2018 pto_search)
    std::vector<Vertex> pto_search(const SearchOrders& so, const std::vector<Vertex>& start, StopRule stop)
    {
        ++res_.pto_searches;
        PseudoTopoOrder t(g_, scc_, so.in_adj, so.out_adj, source_);
        t.incremental = opt_.pto_incremental;
        t.randomize(rng_);
        t.impose(start);
        Weight cur = t.value();
        record(cur);
        Path best_seen = t.best_path(); // with free_first the DP value can decrease
        const bool can_open = opt_.pto_mode == PtoMode::Reimpose || t.has_openable_block(opt_.opener);
        Stopwatch since;
        std::int64_t passes_since = 0;
        std::int64_t sideways_since = 0;
        while (can_open && !deadline_.expired())
        {
            if (stop.seconds_without_improvement > 0 && since.seconds() >= stop.seconds_without_improvement)
                break;
            if (stop.steps_without_improvement > 0 && passes_since >= stop.steps_without_improvement)
                break;
            bool improved = false;
            if (opt_.pto_mode == PtoMode::Open)
                improved = t.open_pass(opt_.opener, rng_, deadline_, res_.open);
            else
            {
                auto p = t.best_path().vertices;
                t.randomize(rng_);
                t.impose(p);
                ++res_.open.passes;
                improved = t.value() > cur;
            }
            if (!improved && res_.open.last_pass_sideways)
            {
                if (++sideways_since >= opt_.sideways_max)
                    break;
                continue; // the path changed: start a new pass without counting a failed one
            }
            ++passes_since;
            if (improved)
            {
                sideways_since = 0;
                cur = t.value();
                record(cur);
                since.reset();
                passes_since = 0;
                if (cur > best_seen.value)
                    best_seen = t.best_path();
            }
            else if (opt_.opener.free_first)
                cur = t.value(); // 2018: the next pass compares with the current (possibly lower) value
        }
        auto bp = t.best_path();
        if (bp.value < best_seen.value)
            bp = best_seen;
        offer(bp.vertices, bp.value);
        return bp.vertices;
    }

    Params next_params(int i)
    {
        if (opt_.params == "2018")
            return params_2018(i, rng_);
        if (opt_.params == "fixed")
            return params_2018(0, rng_);
        if (opt_.params == "random")
            return params_2018(2, rng_);
        throw std::invalid_argument("bad params mode: " + opt_.params);
    }

    // ---- 2018 schedule
    void classic()
    {
        std::vector<Vertex> best;
        Weight bv = kNegInf;
        SearchOrders so;
        auto consider = [&](std::vector<Vertex> p) {
            const Weight v = path_value(p);
            if (!p.empty() && v > bv)
            {
                bv = v;
                best = std::move(p);
            }
        };
        for (int i = 0; i < opt_.param_restarts && !deadline_.expired(); ++i)
        {
            ++res_.restarts;
            so = make_orders(g_, scc_, next_params(i), opt_.ordering, rng_);
            auto f = dfs_forward_full(so);
            if (opt_.pto && opt_.pto_each_dfs && !deadline_.expired())
                f = pto_search(so, f, opt_.pto_stop);
            consider(std::move(f));
            if (opt_.backward && source_ == kNoVertex && !deadline_.expired())
            {
                auto b = dfs_backward_full(so);
                if (!b.empty() && opt_.pto && opt_.pto_each_dfs && !deadline_.expired())
                    b = pto_search(so, b, opt_.pto_stop);
                consider(std::move(b));
            }
        }
        if (opt_.pto && !best.empty() && !deadline_.expired())
            consider(pto_search(so, best, opt_.pto_stop));
    }

    // ---- modern anytime schedule: iterated local search around the PTO neighbourhood
    // starting path from opt_.init_path (original vertex ids); empty if it is not a valid simple path of g_
    std::vector<Vertex> load_init_path() const
    {
        std::ifstream in(opt_.init_path);
        if (!in)
            throw std::invalid_argument("cannot read init_path " + opt_.init_path);
        std::map<long long, Vertex> new_of_old;
        for (Vertex v = 0; v < g_.num_vertices(); ++v)
            new_of_old[static_cast<long long>(old_of_new_[static_cast<std::size_t>(v)])] = v;
        std::vector<Vertex> p;
        long long x = 0;
        while (in >> x)
        {
            auto it = new_of_old.find(x);
            if (it == new_of_old.end())
                return {};
            p.push_back(it->second);
        }
        std::vector<char> seen(static_cast<std::size_t>(g_.num_vertices()), 0);
        for (std::size_t i = 0; i < p.size(); ++i)
        {
            if (seen[static_cast<std::size_t>(p[i])] || (i > 0 && !g_.arc_weight(p[i - 1], p[i])))
                return {};
            seen[static_cast<std::size_t>(p[i])] = 1;
        }
        if (source_ != kNoVertex && (p.empty() || p.front() != source_))
            return {};
        return p;
    }

    // ---- assignment-relaxation opening (ap.hpp): bound, patched path, early stop at the bound
    std::vector<Vertex> ap_opening()
    {
        auto& rep = res_.ap;
        rep.ran = true;
        Stopwatch sw;
        const ApSolution ap = solve_ap(g_, source_, deadline_);
        rep.ap_seconds = sw.seconds();
        rep.solved = ap.solved;
        rep.bids = ap.bids;
        rep.phases = ap.phases;
        if (!ap.solved)
            return {};
        rep.certified = ap.certified;
        rep.bound = ap.bound;
        rep.bound_global = bound_global_;
        rep.path_vertices = ap.path.size();
        rep.cycles = ap.cycles.size();
        for (const auto& c : ap.cycles)
            rep.cycle_vertices += c.size();
        if (ap.certified && bound_global_)
            bound_ = ap.bound;
        Stopwatch sp;
        auto p = patch_ap(g_, ap, source_, &rep.patch);
        rep.patch_seconds = sp.seconds();
        rep.patched = path_value(p);
        offer(p, rep.patched);
        return p;
    }
    [[nodiscard]] bool proven_optimal() const { return best_path_value_ >= bound_; }

    void ils()
    {
        std::vector<Vertex> ap_path;
        if (opt_.ap_init)
        {
            ap_path = ap_opening();
            if (proven_optimal())
                return;
        }
        SearchOrders so = make_orders(g_, scc_, next_params(0), opt_.ordering, rng_);
        std::vector<Vertex> cur = dfs_forward_full(so);
        if (source_ == kNoVertex && opt_.backward && !deadline_.expired())
        {
            auto b = dfs_backward_full(so);
            if (path_value(b) > path_value(cur))
                cur = std::move(b);
        }
        if (!opt_.init_path.empty())
        {
            auto ip = load_init_path();
            if (!ip.empty() && path_value(ip) > path_value(cur))
                cur = std::move(ip);
        }
        if (!ap_path.empty() && path_value(ap_path) > path_value(cur))
            cur = std::move(ap_path);
        offer(cur, path_value(cur));
        if (opt_.pto && !deadline_.expired())
            cur = pto_search(so, cur, opt_.pto_stop);
        Weight cv = path_value(cur);
        StopRule ils_dfs = opt_.dfs_stop;
        if (opt_.ils_dfs_steps > 0)
            ils_dfs = {0, opt_.ils_dfs_steps};
        int round = 1;
        while (!deadline_.expired() && !proven_optimal() && (opt_.ils_rounds == 0 || round <= opt_.ils_rounds))
        {
            ++res_.restarts;
            if (opt_.param_restarts > 1 && round % 4 == 0)
                so = make_orders(g_, scc_, next_params(round), opt_.ordering, rng_);
            ++round;
            // perturbation: keep a random window of the current path
            std::vector<Vertex> p = cur;
            const auto k = p.size();
            if (k >= 4)
            {
                auto keep = std::max<std::size_t>(2, static_cast<std::size_t>(opt_.ils_keep_fraction * static_cast<double>(k)));
                std::size_t off = 0;
                if (source_ == kNoVertex)
                    off = rng_.below(k - keep + 1);
                p.assign(cur.begin() + static_cast<std::ptrdiff_t>(off), cur.begin() + static_cast<std::ptrdiff_t>(off + keep));
            }
            p = dfs_forward(so, p, ils_dfs);
            if (source_ == kNoVertex && !deadline_.expired())
                p = dfs_backward(so, p, ils_dfs);
            if (opt_.pto && !deadline_.expired())
                p = pto_search(so, p, opt_.pto_stop);
            const Weight pv = path_value(p);
            if (pv >= cv)
            {
                cv = pv;
                cur = std::move(p);
            }
        }
    }


    // ---- GAMM: the mutation-based genetic algorithm of Portugal, Antunes and Rocha (2010),
    // re-implemented from the paper's description and adapted to digraphs (see code/README.md).
    struct Individual
    {
        std::vector<Vertex> p;
        Weight value;
    };

    // degree-biased random extension (Portugal et al.: next vertex chosen with probability
    // proportional to its degree); forward uses out-arcs/out-degree, backward in-arcs/in-degree.
    Weight extend_random(std::vector<Vertex>& seq, std::vector<char>& used, bool forward, Vertex forbid)
    {
        Weight gained = 0;
        std::vector<std::pair<Vertex, Weight>> cand;
        std::vector<double> wts;
        while (true)
        {
            const Vertex u = seq.back();
            cand.clear();
            wts.clear();
            double total = 0;
            for (const auto& a : forward ? g_.out(u) : g_.in(u))
            {
                if (used[static_cast<std::size_t>(a.v)] || a.v == forbid)
                    continue;
                const double d = static_cast<double>(forward ? g_.outdegree(a.v) : g_.indegree(a.v));
                cand.emplace_back(a.v, a.w);
                wts.push_back(std::max(d, 0.5));
                total += wts.back();
            }
            forbid = kNoVertex;
            if (cand.empty())
                return gained;
            double r = rng_.uniform01() * total;
            std::size_t k = 0;
            while (k + 1 < cand.size() && r >= wts[k])
                r -= wts[k++];
            seq.push_back(cand[k].first);
            used[static_cast<std::size_t>(cand[k].first)] = 1;
            gained += cand[k].second;
        }
    }

    Individual random_individual(std::vector<char>& used)
    {
        const Vertex s = source_ != kNoVertex ? source_ : random_vertex();
        std::vector<Vertex> f{s};
        used[static_cast<std::size_t>(s)] = 1;
        Weight v = extend_random(f, used, true, kNoVertex);
        if (source_ == kNoVertex)
        {
            std::vector<Vertex> b{s};
            v += extend_random(b, used, false, kNoVertex);
            std::vector<Vertex> p(b.rbegin(), b.rend());
            p.insert(p.end(), f.begin() + 1, f.end());
            f = std::move(p);
        }
        for (Vertex x : f)
            used[static_cast<std::size_t>(x)] = 0;
        return {std::move(f), v};
    }

    // mutation: keep the path up to the cut vertex and regrow from it avoiding the old successor
    bool mutate(const Individual& in, std::size_t pressure, bool at_end, Individual& out, std::vector<char>& used)
    {
        const std::size_t k = in.p.size();
        if (k < 2)
            return false;
        std::size_t back = std::min<std::size_t>(pressure, (k - 1) / 2 + 1);
        std::vector<Vertex> seq;
        Vertex forbid = kNoVertex;
        if (at_end)
        {
            const std::size_t cut = k - 1 - std::min(back, k - 1); // keep p[0..cut]
            seq.assign(in.p.begin(), in.p.begin() + static_cast<std::ptrdiff_t>(cut + 1));
            forbid = in.p[cut + 1];
        }
        else
        {
            const std::size_t cut = std::min(back, k - 1); // keep p[cut..k)
            seq.assign(in.p.rbegin(), in.p.rbegin() + static_cast<std::ptrdiff_t>(k - cut));
            forbid = in.p[cut - 1];
        }
        for (Vertex x : seq)
            used[static_cast<std::size_t>(x)] = 1;
        extend_random(seq, used, at_end, forbid);
        for (Vertex x : seq)
            used[static_cast<std::size_t>(x)] = 0;
        if (!at_end)
            std::ranges::reverse(seq);
        if (seq == in.p)
            return false;
        out.p = std::move(seq);
        out.value = path_value(out.p);
        return true;
    }

    void gamm()
    {
        const std::size_t pop = static_cast<std::size_t>(std::max(4, opt_.ga_population));
        const std::size_t elite = std::max<std::size_t>(1, static_cast<std::size_t>(opt_.ga_elite * static_cast<double>(pop)));
        std::vector<char> used(static_cast<std::size_t>(g_.num_vertices()), 0);
        while (!deadline_.expired())
        {
            ++res_.restarts;
            std::vector<Individual> popu;
            for (std::size_t i = 0; i < pop && !deadline_.expired(); ++i)
            {
                popu.push_back(random_individual(used));
                offer(popu.back().p, popu.back().value);
            }
            Weight best = kNegInf;
            for (const auto& x : popu)
                best = std::max(best, x.value);
            std::size_t pressure = 1;
            int stagnation = 0;
            while (stagnation < opt_.ga_stagnation && !deadline_.expired())
            {
                std::vector<Individual> pool = popu;
                std::size_t successes = 0;
                for (const auto& x : popu)
                {
                    if (deadline_.expired())
                        break;
                    Individual c;
                    for (int side = 0; side < (source_ == kNoVertex ? 2 : 1); ++side)
                        if (mutate(x, pressure, side == 0, c, used))
                        {
                            ++successes;
                            offer(c.p, c.value);
                            pool.push_back(std::move(c));
                        }
                }
                // remove duplicates, elitism, roulette wheel on the rest
                std::ranges::sort(pool, [](const Individual& a, const Individual& b) {
                    return a.value != b.value ? a.value > b.value : a.p < b.p;
                });
                pool.erase(std::unique(pool.begin(), pool.end(), [](const Individual& a, const Individual& b) { return a.p == b.p; }),
                           pool.end());
                std::vector<Individual> next(pool.begin(), pool.begin() + static_cast<std::ptrdiff_t>(std::min(elite, pool.size())));
                std::vector<std::size_t> rest;
                for (std::size_t i = next.size(); i < pool.size(); ++i)
                    rest.push_back(i);
                const Weight minv = pool.empty() ? 0 : pool.back().value;
                while (next.size() < pop && !rest.empty())
                {
                    double total = 0;
                    for (auto i : rest)
                        total += static_cast<double>(pool[i].value - minv + 1);
                    double r = rng_.uniform01() * total;
                    std::size_t k = 0;
                    while (k + 1 < rest.size() && r >= static_cast<double>(pool[rest[k]].value - minv + 1))
                        r -= static_cast<double>(pool[rest[k++]].value - minv + 1);
                    next.push_back(std::move(pool[rest[k]]));
                    rest.erase(rest.begin() + static_cast<std::ptrdiff_t>(k));
                }
                popu = std::move(next);
                const Weight nb = popu.empty() ? kNegInf : popu.front().value;
                if (nb > best)
                {
                    best = nb;
                    stagnation = 0;
                }
                else
                {
                    ++stagnation;
                    ++pressure;
                }
                if (successes == 0)
                    pressure = 1;
            }
        }
    }

    // ---- baselines
    void rdfs()
    {
        OrderedAdjacency adj = plain_out_adjacency(g_);
        const std::int64_t steps = opt_.rdfs_steps > 0 ? opt_.rdfs_steps : 10LL * g_.num_vertices();
        InPlaceDfs d(g_.num_vertices(), adj);
        while (!deadline_.expired())
        {
            ++res_.restarts;
            for (Vertex v = 0; v < g_.num_vertices(); ++v)
            {
                auto first = adj.arcs.begin() + static_cast<std::ptrdiff_t>(adj.off[static_cast<std::size_t>(v)]);
                auto last = adj.arcs.begin() + static_cast<std::ptrdiff_t>(adj.off[static_cast<std::size_t>(v) + 1]);
                rng_.shuffle(first, last);
            }
            const Vertex s = source_ != kNoVertex ? source_ : random_vertex();
            d.reset({s});
            d.greedy_extend();
            auto p = d.run({0, steps}, deadline_, res_.dfs, sink());
            offer(p, path_value(p));
        }
    }

    void greedy()
    {
        const Vertex n = g_.num_vertices();
        std::vector<char> used(static_cast<std::size_t>(n), 0);
        std::vector<Vertex> fwd;
        std::vector<Vertex> bwd;
        const bool unit = g_.is_unit_weight();
        auto free_out = [&](Vertex v) {
            int c = 0;
            for (const auto& a : g_.out(v))
                c += used[static_cast<std::size_t>(a.v)] ? 0 : 1;
            return c;
        };
        auto free_in = [&](Vertex v) {
            int c = 0;
            for (const auto& a : g_.in(v))
                c += used[static_cast<std::size_t>(a.v)] ? 0 : 1;
            return c;
        };
        while (!deadline_.expired())
        {
            ++res_.restarts;
            const Vertex s = source_ != kNoVertex ? source_ : random_vertex();
            fwd.assign(1, s);
            bwd.clear();
            used[static_cast<std::size_t>(s)] = 1;
            Weight val = 0;
            // extend the tail, then the head; Warnsdorff rule (fewest free continuations),
            // heaviest arc first on weighted graphs; ties uniformly at random
            for (int dir = 0; dir < (source_ == kNoVertex ? 2 : 1); ++dir)
                while (true)
                {
                    const Vertex end = dir == 0 ? fwd.back() : (bwd.empty() ? fwd.front() : bwd.back());
                    auto arcs = dir == 0 ? g_.out(end) : g_.in(end);
                    Vertex pick = kNoVertex;
                    Weight pw = 0;
                    long long best_key = 0;
                    std::uint64_t ties = 0;
                    for (const auto& a : arcs)
                    {
                        if (used[static_cast<std::size_t>(a.v)])
                            continue;
                        // dead ends (no free continuation) are taken only as a last resort
                        long long deg = dir == 0 ? free_out(a.v) : free_in(a.v);
                        if (deg == 0)
                            deg = 1LL << 40;
                        const long long key = unit ? deg : -a.w * (1LL << 42) + deg;
                        if (pick == kNoVertex || key < best_key)
                        {
                            pick = a.v;
                            pw = a.w;
                            best_key = key;
                            ties = 1;
                        }
                        else if (key == best_key && rng_.below(++ties) == 0)
                        {
                            pick = a.v;
                            pw = a.w;
                        }
                    }
                    if (pick == kNoVertex)
                        break;
                    used[static_cast<std::size_t>(pick)] = 1;
                    val += pw;
                    (dir == 0 ? fwd : bwd).push_back(pick);
                }
            std::vector<Vertex> p(bwd.rbegin(), bwd.rend());
            p.insert(p.end(), fwd.begin(), fwd.end());
            for (Vertex v : p)
                used[static_cast<std::size_t>(v)] = 0;
            offer(p, val);
        }
    }

    const SolverOptions& opt_;
    Rng rng_;
    Stopwatch clock_;
    Deadline deadline_;
    Digraph g_;
    std::vector<Vertex> old_of_new_;
    Vertex source_ = kNoVertex;
    SccInfo scc_;
    SolveResult res_;
    std::vector<Vertex> best_;
    Weight best_value_ = kNegInf;
    Weight best_path_value_ = kNegInf;
    bool bound_global_ = false;
    Weight bound_ = std::numeric_limits<Weight>::max(); // certified upper bound of the whole problem
};

} // namespace

SolveResult solve(const Digraph& g, const SolverOptions& opt)
{
    Solver s(g, opt);
    return s.run();
}

} // namespace lsp
