#include "lsp/structure.hpp"

#include <numeric>

namespace lsp {

SccInfo compute_scc(const Digraph& g)
{
    const Vertex n = g.num_vertices();
    const auto un = static_cast<std::size_t>(n);
    std::vector<std::int32_t> index(un, -1);
    std::vector<std::int32_t> low(un, 0);
    std::vector<char> on_stack(un, 0);
    std::vector<Vertex> stack;
    std::vector<std::int32_t> comp(un, -1);
    std::int32_t next_index = 0;
    std::int32_t num_comp = 0;

    struct Frame
    {
        Vertex v;
        std::size_t arc;
    };
    std::vector<Frame> call;
    for (Vertex root = 0; root < n; ++root)
    {
        if (index[static_cast<std::size_t>(root)] != -1)
            continue;
        call.push_back({root, 0});
        index[static_cast<std::size_t>(root)] = low[static_cast<std::size_t>(root)] = next_index++;
        stack.push_back(root);
        on_stack[static_cast<std::size_t>(root)] = 1;
        while (!call.empty())
        {
            Frame& f = call.back();
            const auto out = g.out(f.v);
            if (f.arc < out.size())
            {
                const Vertex w = out[f.arc++].v;
                const auto uw = static_cast<std::size_t>(w);
                if (index[uw] == -1)
                {
                    index[uw] = low[uw] = next_index++;
                    stack.push_back(w);
                    on_stack[uw] = 1;
                    call.push_back({w, 0});
                }
                else if (on_stack[uw])
                {
                    low[static_cast<std::size_t>(f.v)] = std::min(low[static_cast<std::size_t>(f.v)], index[uw]);
                }
                continue;
            }
            const Vertex v = f.v;
            const auto uv = static_cast<std::size_t>(v);
            call.pop_back();
            if (!call.empty())
            {
                const auto up = static_cast<std::size_t>(call.back().v);
                low[up] = std::min(low[up], low[uv]);
            }
            if (low[uv] == index[uv])
            {
                while (true)
                {
                    const Vertex w = stack.back();
                    stack.pop_back();
                    on_stack[static_cast<std::size_t>(w)] = 0;
                    comp[static_cast<std::size_t>(w)] = num_comp;
                    if (w == v)
                        break;
                }
                ++num_comp;
            }
        }
    }
    // Tarjan emits SCCs in reverse topological order: flip ids.
    SccInfo info;
    info.num_components = num_comp;
    info.comp.resize(un);
    for (std::size_t v = 0; v < un; ++v)
        info.comp[v] = num_comp - 1 - comp[v];

    const auto uc = static_cast<std::size_t>(num_comp);
    info.start.assign(uc + 1, 0);
    for (std::size_t v = 0; v < un; ++v)
        ++info.start[static_cast<std::size_t>(info.comp[v]) + 1];
    for (std::size_t c = 0; c < uc; ++c)
        info.start[c + 1] += info.start[c];
    info.block.resize(un);
    std::vector<std::int32_t> fill(info.start.begin(), info.start.end() - 1);
    for (Vertex v = 0; v < n; ++v)
        info.block[static_cast<std::size_t>(fill[static_cast<std::size_t>(info.comp[static_cast<std::size_t>(v)])]++)] = v;

    info.rank_out.assign(uc, 0);
    info.rank_in.assign(uc, 0);
    // rank_in: process components in topological order
    for (std::size_t c = 0; c < uc; ++c)
        for (auto i = info.start[c]; i < info.start[c + 1]; ++i)
        {
            const Vertex v = info.block[static_cast<std::size_t>(i)];
            for (const auto& a : g.in(v))
            {
                const auto cu = static_cast<std::size_t>(info.comp[static_cast<std::size_t>(a.v)]);
                if (cu != c)
                    info.rank_in[c] = std::max(info.rank_in[c], info.rank_in[cu] + 1);
            }
        }
    for (std::size_t cc = uc; cc-- > 0;)
        for (auto i = info.start[cc]; i < info.start[cc + 1]; ++i)
        {
            const Vertex v = info.block[static_cast<std::size_t>(i)];
            for (const auto& a : g.out(v))
            {
                const auto cw = static_cast<std::size_t>(info.comp[static_cast<std::size_t>(a.v)]);
                if (cw != cc)
                    info.rank_out[cc] = std::max(info.rank_out[cc], info.rank_out[cw] + 1);
            }
        }
    return info;
}

std::vector<std::int32_t> weak_components(const Digraph& g, std::int32_t& count)
{
    const auto un = static_cast<std::size_t>(g.num_vertices());
    std::vector<std::int32_t> parent(un);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](std::int32_t x) {
        while (parent[static_cast<std::size_t>(x)] != x)
        {
            parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
            x = parent[static_cast<std::size_t>(x)];
        }
        return x;
    };
    for (Vertex u = 0; u < g.num_vertices(); ++u)
        for (const auto& a : g.out(u))
        {
            const auto ru = find(u);
            const auto rv = find(a.v);
            if (ru != rv)
                parent[static_cast<std::size_t>(std::max(ru, rv))] = std::min(ru, rv);
        }
    std::vector<std::int32_t> id(un, -1);
    std::vector<std::int32_t> wcc(un);
    count = 0;
    for (std::size_t v = 0; v < un; ++v)
    {
        const auto r = static_cast<std::size_t>(find(static_cast<std::int32_t>(v)));
        if (id[r] == -1)
            id[r] = count++;
        wcc[v] = id[r];
    }
    return wcc;
}

std::vector<char> reachable_from(const Digraph& g, Vertex s)
{
    std::vector<char> seen(static_cast<std::size_t>(g.num_vertices()), 0);
    std::vector<Vertex> st{s};
    seen[static_cast<std::size_t>(s)] = 1;
    while (!st.empty())
    {
        const Vertex u = st.back();
        st.pop_back();
        for (const auto& a : g.out(u))
            if (!seen[static_cast<std::size_t>(a.v)])
            {
                seen[static_cast<std::size_t>(a.v)] = 1;
                st.push_back(a.v);
            }
    }
    return seen;
}

} // namespace lsp
