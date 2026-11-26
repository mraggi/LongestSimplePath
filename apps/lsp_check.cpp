// Verifies a path file against a graph: simple, every step is an arc, optional claimed value.
//   lsp_check GRAPH PATHFILE [--format F] [--value V] [--source S]
#include "lsp/io.hpp"

#include <iostream>

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        std::cerr << "usage: lsp_check GRAPH PATHFILE [--format F] [--value V] [--source S]\n";
        return 2;
    }
    lsp::GraphFormat fmt = lsp::GraphFormat::Auto;
    std::optional<lsp::Weight> claimed;
    std::optional<lsp::Vertex> source;
    for (int i = 3; i + 1 < argc; i += 2)
    {
        const std::string a = argv[i];
        if (a == "--format")
            fmt = lsp::parse_format(argv[i + 1]);
        else if (a == "--value")
            claimed = std::stoll(argv[i + 1]);
        else if (a == "--source")
            source = static_cast<lsp::Vertex>(std::stoll(argv[i + 1]));
    }
    try
    {
        const auto g = lsp::read_graph(argv[1], fmt);
        const auto p = lsp::read_path(argv[2]);
        const auto r = lsp::check_path(g, p, claimed, source);
        std::cout << (r.ok ? "ok" : "FAIL: " + r.message) << " value=" << r.recomputed << " vertices=" << p.size() << '\n';
        return r.ok ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "error: " << e.what() << '\n';
        return 2;
    }
}
