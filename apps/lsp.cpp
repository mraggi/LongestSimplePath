// Command line interface.
//
//   lsp solve GRAPH [--format F] [--method M] [--time T] [--seed S] [--source V]
//                   [--set key=value]... [--path-out FILE] [--print-path]
//   lsp info  GRAPH [--format F]
//   lsp convert GRAPH OUT [--format F]      (writes the 0-based edge-list format)
//
// `solve` prints one JSON object on stdout (see code/README.md).
#include "lsp/io.hpp"
#include "lsp/solver.hpp"

#include <fstream>
#include <iostream>

namespace {

void usage()
{
    std::cerr << "usage:\n"
                 "  lsp solve GRAPH [--format F] [--preset 2018|new|newap] [--method classic|ils|rdfs|greedy|gamm] [--time T] [--seed S]\n"
                 "            [--source V] [--set key=value]... [--path-out FILE] [--print-path]\n"
                 "  lsp info GRAPH [--format F]\n"
                 "  lsp ap GRAPH [--format F] [--source V] [--path-out FILE]   (assignment bound + patched path only)\n"
                 "  lsp convert GRAPH OUT [--format F]\n";
}

std::string json_escape(const std::string& s)
{
    std::string r;
    for (char c : s)
    {
        if (c == '"' || c == '\\')
            r += '\\';
        r += c;
    }
    return r;
}

int cmd_info(const lsp::Digraph& g)
{
    auto scc = lsp::compute_scc(g);
    std::int32_t nw = 0;
    auto w = lsp::weak_components(g, nw);
    std::int32_t biggest = 0;
    std::int32_t nontrivial = 0;
    for (std::int32_t c = 0; c < scc.num_components; ++c)
    {
        biggest = std::max(biggest, scc.size(c));
        nontrivial += scc.size(c) > 1 ? 1 : 0;
    }
    std::int32_t maxrank = 0;
    for (auto r : scc.rank_out)
        maxrank = std::max(maxrank, r);
    std::cout << "{\"n\": " << g.num_vertices() << ", \"m\": " << g.num_edges()
              << ", \"unit_weight\": " << (g.is_unit_weight() ? "true" : "false") << ", \"wcc\": " << nw
              << ", \"scc\": " << scc.num_components << ", \"nontrivial_scc\": " << nontrivial
              << ", \"largest_scc\": " << biggest << ", \"condensation_depth\": " << maxrank << "}\n";
    return 0;
}

int cmd_ap(const lsp::Digraph& input, lsp::Vertex source, const std::string& path_out)
{
    // same preprocessing as the rooted solver: vertices reachable from the source, no arcs into it
    lsp::Digraph g = input;
    std::vector<lsp::Vertex> old_of_new;
    lsp::Vertex s = lsp::kNoVertex;
    if (source != lsp::kNoVertex)
    {
        g = input.induced(lsp::reachable_from(input, source), old_of_new);
        for (lsp::Vertex v = 0; v < g.num_vertices(); ++v)
            if (old_of_new[static_cast<std::size_t>(v)] == source)
                s = v;
        g = g.without_in_arcs(s);
    }
    lsp::Stopwatch c1;
    const auto ap = lsp::solve_ap(g, s);
    const double t1 = c1.seconds();
    lsp::Stopwatch c2;
    lsp::PatchStats ps;
    auto p = lsp::patch_ap(g, ap, s, &ps);
    const double t2 = c2.seconds();
    if (!old_of_new.empty())
        for (auto& v : p)
            v = old_of_new[static_cast<std::size_t>(v)];
    const auto chk = lsp::check_path(input, p, std::nullopt,
                                     source != lsp::kNoVertex ? std::optional<lsp::Vertex>(source) : std::nullopt);
    std::size_t cv = 0;
    for (const auto& c : ap.cycles)
        cv += c.size();
    std::cout << "{\"n\": " << input.num_vertices() << ", \"m\": " << input.num_edges() << ", \"ap_bound\": " << ap.bound
              << ", \"ap_certified\": " << (ap.certified ? "true" : "false") << ", \"ap_path_vertices\": " << ap.path.size()
              << ", \"ap_cycles\": " << ap.cycles.size() << ", \"ap_cycle_vertices\": " << cv << ", \"patched\": " << chk.recomputed
              << ", \"check\": \"" << (chk.ok ? std::string("ok") : json_escape(chk.message)) << "\", \"ap_seconds\": " << t1
              << ", \"patch_seconds\": " << t2 << ", \"ap_bids\": " << ap.bids << ", \"ap_phases\": " << ap.phases
              << ", \"patch_merged\": " << ps.merged << ", \"patch_lossy\": " << ps.lossy << ", \"patch_left\": " << ps.left
              << ", \"patch_extended\": " << ps.extended << "}\n";
    if (!path_out.empty())
    {
        std::ofstream out(path_out);
        lsp::write_path(out, p);
    }
    return chk.ok ? 0 : 1;
}

} // namespace

int main(int argc, char** argv)
{
    std::ios_base::sync_with_stdio(false);
    if (argc < 3)
    {
        usage();
        return 2;
    }
    const std::string cmd = argv[1];
    const std::string file = argv[2];
    lsp::GraphFormat fmt = lsp::GraphFormat::Auto;
    lsp::SolverOptions opt;
    std::string path_out;
    std::string convert_out;
    bool print_path = false;
    int i = 3;
    if (cmd == "convert")
    {
        if (argc < 4)
        {
            usage();
            return 2;
        }
        convert_out = argv[3];
        i = 4;
    }
    try
    {
        for (; i < argc; ++i)
        {
            const std::string a = argv[i];
            auto next = [&]() -> std::string {
                if (i + 1 >= argc)
                    throw std::invalid_argument("missing value after " + a);
                return argv[++i];
            };
            if (a == "--format")
                fmt = lsp::parse_format(next());
            else if (a == "--preset")
                opt.set("preset", next());
            else if (a == "--method")
                opt.set("method", next());
            else if (a == "--time")
                opt.set("time", next());
            else if (a == "--seed")
                opt.set("seed", next());
            else if (a == "--source")
                opt.set("source", next());
            else if (a == "--set")
            {
                const std::string kv = next();
                const auto eq = kv.find('=');
                if (eq == std::string::npos)
                    throw std::invalid_argument("--set expects key=value");
                opt.set(kv.substr(0, eq), kv.substr(eq + 1));
            }
            else if (a == "--path-out")
                path_out = next();
            else if (a == "--print-path")
                print_path = true;
            else
                throw std::invalid_argument("unknown argument " + a);
        }
        lsp::Stopwatch read_clock;
        const lsp::Digraph g = lsp::read_graph(file, fmt);
        const double read_seconds = read_clock.seconds();
        if (cmd == "info")
            return cmd_info(g);
        if (cmd == "ap")
            return cmd_ap(g, opt.source, path_out);
        if (cmd == "convert")
        {
            std::ofstream out(convert_out);
            lsp::write_edgelist(out, g);
            return 0;
        }
        if (cmd != "solve")
        {
            usage();
            return 2;
        }
        const auto res = lsp::solve(g, opt);
        const auto chk = lsp::check_path(g, res.best.vertices, res.best.value,
                                         opt.source != lsp::kNoVertex ? std::optional<lsp::Vertex>(opt.source) : std::nullopt);
        if (!path_out.empty())
        {
            std::ofstream out(path_out);
            lsp::write_path(out, res.best.vertices);
        }
        auto& o = std::cout;
        o << "{\"graph\": \"" << json_escape(file) << "\", \"n\": " << g.num_vertices() << ", \"m\": " << g.num_edges()
          << ", \"work_n\": " << res.work_vertices << ", \"work_m\": " << res.work_edges
          << ", \"value\": " << res.best.value << ", \"num_vertices_path\": " << res.best.vertices.size()
          << ", \"check\": \"" << (chk.ok ? std::string("ok") : json_escape(chk.message)) << "\""
          << ", \"read_seconds\": " << read_seconds << ", \"preprocess_seconds\": " << res.preprocess_seconds
          << ", \"total_seconds\": " << res.total_seconds << ", \"dfs_steps\": " << res.dfs.steps
          << ", \"openings\": " << res.open.openings << ", \"open_passes\": " << res.open.passes
          << ", \"open_improvements\": " << res.open.improvements << ", \"dp_vertices\": " << res.open.dp_vertices
          << ", \"sort_ops\": " << res.open.sort_ops << ", \"pto_searches\": " << res.pto_searches
          << ", \"restarts\": " << res.restarts;
        if (res.ap.ran)
        {
            const auto& a = res.ap;
            o << ", \"ap_solved\": " << (a.solved ? "true" : "false") << ", \"ap_certified\": " << (a.certified ? "true" : "false")
              << ", \"ap_bound_global\": " << (a.bound_global ? "true" : "false") << ", \"ap_bound\": " << a.bound
              << ", \"ap_path_vertices\": " << a.path_vertices << ", \"ap_cycles\": " << a.cycles
              << ", \"ap_cycle_vertices\": " << a.cycle_vertices << ", \"ap_patched\": " << a.patched
              << ", \"ap_seconds\": " << a.ap_seconds << ", \"patch_seconds\": " << a.patch_seconds << ", \"ap_bids\": " << a.bids
              << ", \"ap_phases\": " << a.phases << ", \"patch_merged\": " << a.patch.merged << ", \"patch_lossy\": " << a.patch.lossy
              << ", \"patch_left\": " << a.patch.left << ", \"patch_extended\": " << a.patch.extended;
        }
        o << ", \"proven_optimal\": " << (res.proven_optimal ? "true" : "false") << ", \"options\": {";
        bool first = true;
        for (const auto& [k, v] : opt.describe())
        {
            o << (first ? "" : ", ") << '"' << k << "\": \"" << json_escape(v) << '"';
            first = false;
        }
        o << "}, \"trajectory\": [";
        for (std::size_t t = 0; t < res.trajectory.size(); ++t)
            o << (t ? ", " : "") << '[' << res.trajectory[t].time << ", " << res.trajectory[t].value << ']';
        o << ']';
        if (print_path)
        {
            o << ", \"path\": [";
            for (std::size_t t = 0; t < res.best.vertices.size(); ++t)
                o << (t ? ", " : "") << res.best.vertices[t];
            o << ']';
        }
        o << "}\n";
        return chk.ok ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "error: " << e.what() << '\n';
        return 2;
    }
}
