#include "lsp/io.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace lsp {

GraphFormat parse_format(const std::string& name)
{
    if (name == "auto")
        return GraphFormat::Auto;
    if (name == "edgelist")
        return GraphFormat::EdgeList;
    if (name == "adjlist")
        return GraphFormat::AdjList;
    if (name == "snap")
        return GraphFormat::Snap;
    if (name == "dimacs")
        return GraphFormat::Dimacs;
    throw std::invalid_argument("unknown graph format: " + name);
}

namespace {

bool ends_with(const std::string& s, const std::string& suf)
{
    return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

Digraph read_edgelist(std::istream& in)
{
    std::string line;
    long long n = -1;
    long long m = -1;
    std::vector<Edge> edges;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#' || line[0] == '%')
            continue;
        std::istringstream ss(line);
        if (n < 0)
        {
            ss >> n >> m;
            if (!ss)
                throw std::runtime_error("edgelist: bad header");
            edges.reserve(static_cast<std::size_t>(m));
            continue;
        }
        long long u = 0;
        long long v = 0;
        long long w = 1;
        ss >> u >> v;
        if (!ss)
            throw std::runtime_error("edgelist: bad edge line: " + line);
        if (!(ss >> w))
            w = 1;
        edges.push_back({static_cast<Vertex>(u), static_cast<Vertex>(v), static_cast<Weight>(w)});
    }
    if (n < 0)
        throw std::runtime_error("edgelist: empty input");
    if (static_cast<long long>(edges.size()) != m)
        throw std::runtime_error("edgelist: header says " + std::to_string(m) + " edges, read " +
                                 std::to_string(edges.size()));
    return {static_cast<Vertex>(n), std::move(edges)};
}

Digraph read_adjlist(std::istream& in)
{
    std::string line;
    std::vector<std::vector<long long>> rows;
    while (std::getline(in, line))
    {
        std::istringstream ss(line);
        std::vector<long long> row;
        long long x = 0;
        while (ss >> x)
            row.push_back(x);
        rows.push_back(std::move(row));
    }
    while (!rows.empty() && rows.back().empty())
        rows.pop_back();
    const auto n = static_cast<Vertex>(rows.size());
    std::vector<Edge> edges;
    for (Vertex u = 0; u < n; ++u)
        for (long long x : rows[static_cast<std::size_t>(u)])
        {
            if (x < 1 || x > n)
                throw std::runtime_error("adjlist: neighbour out of range");
            edges.push_back({u, static_cast<Vertex>(x - 1), 1});
        }
    return {n, std::move(edges)};
}

Digraph read_snap(std::istream& in)
{
    std::string line;
    std::unordered_map<long long, Vertex> id;
    std::vector<std::string> labels;
    std::vector<Edge> edges;
    auto get = [&](long long x) {
        auto [it, inserted] = id.try_emplace(x, static_cast<Vertex>(labels.size()));
        if (inserted)
            labels.push_back(std::to_string(x));
        return it->second;
    };
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#' || line[0] == '%')
            continue;
        std::istringstream ss(line);
        long long u = 0;
        long long v = 0;
        long long w = 1;
        ss >> u >> v;
        if (!ss)
            continue;
        if (!(ss >> w))
            w = 1;
        const Vertex a = get(u);
        const Vertex b = get(v);
        edges.push_back({a, b, static_cast<Weight>(w)});
    }
    Digraph g(static_cast<Vertex>(labels.size()), std::move(edges));
    g.labels = std::move(labels);
    return g;
}

Digraph read_dimacs(std::istream& in)
{
    std::string line;
    long long n = -1;
    std::vector<Edge> edges;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == 'c')
            continue;
        std::istringstream ss(line);
        char tag = 0;
        ss >> tag;
        if (tag == 'p')
        {
            std::string kind;
            long long m = 0;
            ss >> kind >> n >> m;
            edges.reserve(static_cast<std::size_t>(m));
        }
        else if (tag == 'a')
        {
            long long u = 0;
            long long v = 0;
            long long w = 1;
            ss >> u >> v >> w;
            edges.push_back({static_cast<Vertex>(u - 1), static_cast<Vertex>(v - 1), static_cast<Weight>(w)});
        }
    }
    if (n < 0)
        throw std::runtime_error("dimacs: missing p line");
    return {static_cast<Vertex>(n), std::move(edges)};
}

} // namespace

Digraph read_graph(std::istream& in, GraphFormat format)
{
    switch (format)
    {
    case GraphFormat::EdgeList:
        return read_edgelist(in);
    case GraphFormat::AdjList:
        return read_adjlist(in);
    case GraphFormat::Snap:
        return read_snap(in);
    case GraphFormat::Dimacs:
        return read_dimacs(in);
    case GraphFormat::Auto:
        break;
    }
    throw std::invalid_argument("read_graph(istream): format must be explicit");
}

Digraph read_graph(const std::string& filename, GraphFormat format)
{
    std::ifstream in(filename);
    if (!in)
        throw std::runtime_error("cannot open " + filename);
    if (format == GraphFormat::Auto)
    {
        if (ends_with(filename, ".adj") || filename.find("GoGame") != std::string::npos)
            format = GraphFormat::AdjList;
        else if (ends_with(filename, ".gr") || ends_with(filename, ".dimacs"))
            format = GraphFormat::Dimacs;
        else if (ends_with(filename, ".snap"))
            format = GraphFormat::Snap;
        else
            format = GraphFormat::EdgeList;
    }
    return read_graph(in, format);
}

void write_edgelist(std::ostream& out, const Digraph& g)
{
    out << g.num_vertices() << ' ' << g.num_edges() << '\n';
    for (Vertex u = 0; u < g.num_vertices(); ++u)
        for (const auto& a : g.out(u))
        {
            out << u << ' ' << a.v;
            if (!g.is_unit_weight())
                out << ' ' << a.w;
            out << '\n';
        }
}

std::vector<Vertex> read_path(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in)
        throw std::runtime_error("cannot open " + filename);
    std::vector<Vertex> p;
    long long x = 0;
    while (in >> x)
        p.push_back(static_cast<Vertex>(x));
    return p;
}

void write_path(std::ostream& out, const std::vector<Vertex>& p)
{
    for (std::size_t i = 0; i < p.size(); ++i)
        out << p[i] << (i + 1 == p.size() ? '\n' : ' ');
}

} // namespace lsp
