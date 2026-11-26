#include "lsp/pto.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>

namespace lsp {

PseudoTopoOrder::PseudoTopoOrder(const Digraph& g, const SccInfo& scc, const OrderedAdjacency& in_adj,
                                 const OrderedAdjacency& out_adj, Vertex source)
    : g_(g), scc_(scc), in_adj_(in_adj), out_adj_(out_adj), source_(source), n_(g.num_vertices())
{
    const auto un = static_cast<std::size_t>(n_);
    order_ = scc.block;
    pos_.resize(un);
    for (std::size_t i = 0; i < un; ++i)
        pos_[static_cast<std::size_t>(order_[i])] = static_cast<std::int32_t>(i);
    dp_.assign(un, 0);
    parent_.assign(un, -1);
    prefbest_.assign(un, 0);
    dirty_ = 0;
    if (source_ != kNoVertex && g.indegree(source_) != 0)
        throw std::invalid_argument("rooted PTO: the source must have no in-arcs");
}

void PseudoTopoOrder::randomize(Rng& rng)
{
    for (std::int32_t c = 0; c < scc_.num_components; ++c)
    {
        const auto a = scc_.start[static_cast<std::size_t>(c)];
        const auto d = scc_.start[static_cast<std::size_t>(c) + 1];
        if (d - a < 2)
            continue;
        rng.shuffle(order_.begin() + a, order_.begin() + d);
        for (auto i = a; i < d; ++i)
            pos_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])] = i;
    }
    dirty_ = 0;
}

void PseudoTopoOrder::impose(const std::vector<Vertex>& path)
{
    scratch_.clear();
    for (Vertex v : path)
        scratch_.push_back(pos_[static_cast<std::size_t>(v)]);
    std::ranges::sort(scratch_);
    for (std::size_t k = 0; k < path.size(); ++k)
    {
        const auto i = scratch_[k];
        if (order_[static_cast<std::size_t>(i)] != path[k])
        {
            order_[static_cast<std::size_t>(i)] = path[k];
            pos_[static_cast<std::size_t>(path[k])] = i;
            touch(i);
        }
    }
}

void PseudoTopoOrder::swap_positions(std::int32_t i, std::int32_t j)
{
    std::swap(order_[static_cast<std::size_t>(i)], order_[static_cast<std::size_t>(j)]);
    pos_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])] = i;
    pos_[static_cast<std::size_t>(order_[static_cast<std::size_t>(j)])] = j;
    touch(std::min(i, j));
}

void PseudoTopoOrder::evaluate(OpenStats* st)
{
    if (!incremental)
        dirty_ = 0;
    for (std::int32_t i = dirty_; i < n_; ++i)
    {
        const auto ui = static_cast<std::size_t>(i);
        const Vertex u = order_[ui];
        Weight best = (source_ == kNoVertex || u == source_) ? 0 : kNegInf;
        std::int32_t par = -1;
        for (const Arc& a : in_adj_[u])
        {
            const auto j = pos_[static_cast<std::size_t>(a.v)];
            if (j >= i)
                continue;
            const Weight dj = dp_[static_cast<std::size_t>(j)];
            if (dj == kNegInf)
                continue;
            const Weight cand = dj + a.w;
            if (cand > best)
            {
                best = cand;
                par = j;
            }
        }
        dp_[ui] = best;
        parent_[ui] = par;
        if (i == 0 || best > dp_[static_cast<std::size_t>(prefbest_[ui - 1])])
            prefbest_[ui] = i;
        else
            prefbest_[ui] = prefbest_[ui - 1];
    }
    if (st)
        st->dp_vertices += n_ - dirty_;
    dirty_ = n_;
}

Weight PseudoTopoOrder::value()
{
    if (n_ == 0)
        return 0;
    evaluate();
    const Weight v = dp_[static_cast<std::size_t>(prefbest_[static_cast<std::size_t>(n_ - 1)])];
    return v == kNegInf ? 0 : v;
}

std::vector<std::int32_t> PseudoTopoOrder::path_positions()
{
    evaluate();
    std::vector<std::int32_t> ps;
    if (n_ == 0)
        return ps;
    std::int32_t i = prefbest_[static_cast<std::size_t>(n_ - 1)];
    while (i != -1)
    {
        ps.push_back(i);
        i = parent_[static_cast<std::size_t>(i)];
    }
    std::ranges::reverse(ps);
    return ps;
}

Path PseudoTopoOrder::best_path()
{
    Path p;
    for (auto i : path_positions())
        p.vertices.push_back(order_[static_cast<std::size_t>(i)]);
    p.value = value();
    return p;
}

bool PseudoTopoOrder::is_strong_pto() const
{
    for (std::int32_t i = 0; i < n_; ++i)
    {
        const Vertex v = order_[static_cast<std::size_t>(i)];
        const auto c = scc_.comp[static_cast<std::size_t>(v)];
        if (i < scc_.start[static_cast<std::size_t>(c)] || i >= scc_.start[static_cast<std::size_t>(c) + 1])
            return false;
        if (pos_[static_cast<std::size_t>(v)] != i)
            return false;
    }
    return true;
}

bool PseudoTopoOrder::has_openable_block(const OpenerOptions& opt) const
{
    for (std::int32_t c = 0; c < scc_.num_components; ++c)
        if (scc_.size(c) > opt.skip_block_size && scc_.size(c) >= opt.min_free && scc_.size(c) >= 2)
            return true;
    return false;
}

// 2018 PseudoTopoOrder::transfer: the block [a,d) holds h0 = b-a path vertices, then the free
// vertices [b,c), then the remaining path vertices [c,d). Moves path vertices across the gap
// with transpositions until exactly h of them precede it (this rotates the free vertices).
void PseudoTopoOrder::transfer(std::int32_t a, std::int32_t b, std::int32_t c, std::int32_t h)
{
    while (h > b - a)
    {
        swap_positions(b, c);
        ++b;
        ++c;
    }
    while (h < b - a)
    {
        --b;
        --c;
        swap_positions(b, c);
    }
}

std::int32_t PseudoTopoOrder::earlier_out_neighbour(std::int32_t lo, std::int32_t r, BOrderRule rule, Rng& rng)
{
    const Vertex u = order_[static_cast<std::size_t>(r)];
    if (rule == BOrderRule::HeuristicFirst)
    {
        for (const Arc& a : out_adj_[u])
        {
            const auto j = pos_[static_cast<std::size_t>(a.v)];
            if (lo <= j && j < r)
                return j;
        }
        return -1;
    }
    // uniform among the candidates (reservoir sampling)
    std::int32_t chosen = -1;
    std::uint64_t seen = 0;
    for (const Arc& a : out_adj_[u])
    {
        const auto j = pos_[static_cast<std::size_t>(a.v)];
        if (lo <= j && j < r)
        {
            ++seen;
            if (rng.below(seen) == 0)
                chosen = j;
        }
    }
    return chosen;
}

void PseudoTopoOrder::order_free_dfs(std::int32_t b, std::int32_t c, Vertex anchor, Rng& rng)
{
    // stamps: 2*cur+0 = in B and unvisited, 2*cur+1 = visited
    if (stamp_.size() != static_cast<std::size_t>(n_))
        stamp_.assign(static_cast<std::size_t>(n_), 0);
    if (cur_stamp_ > 0x7ffffff0U)
    {
        std::ranges::fill(stamp_, 0);
        cur_stamp_ = 0;
    }
    ++cur_stamp_;
    const std::uint32_t in_b = 2 * cur_stamp_;
    const std::uint32_t seen = 2 * cur_stamp_ + 1;
    for (auto i = b; i < c; ++i)
        stamp_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])] = in_b;
    roots_.clear();
    // preferred roots: free out-neighbours of the path vertex right before the gap
    if (anchor != kNoVertex)
    {
        for (const Arc& a : out_adj_[anchor])
            if (stamp_[static_cast<std::size_t>(a.v)] == in_b)
                roots_.push_back(a.v);
        rng.shuffle(roots_.begin(), roots_.end());
    }
    const std::size_t pref = roots_.size();
    for (auto i = b; i < c; ++i)
        roots_.push_back(order_[static_cast<std::size_t>(i)]);
    rng.shuffle(roots_.begin() + static_cast<std::ptrdiff_t>(pref), roots_.end());
    post_.clear();
    for (Vertex r : roots_)
    {
        if (stamp_[static_cast<std::size_t>(r)] != in_b)
            continue;
        stamp_[static_cast<std::size_t>(r)] = seen;
        const auto deg = static_cast<std::uint32_t>(out_adj_[r].size());
        dfs_stack_.push_back({r, 0, deg ? static_cast<std::uint32_t>(rng.below(deg)) : 0});
        while (!dfs_stack_.empty())
        {
            DfsFrame& f = dfs_stack_.back();
            const auto nb = out_adj_[f.v];
            bool pushed = false;
            while (f.k < nb.size())
            {
                std::size_t idx = f.k + f.off; // both < size: one conditional subtraction instead of a division
                if (idx >= nb.size())
                    idx -= nb.size();
                const Vertex w = nb[idx].v;
                ++f.k;
                if (stamp_[static_cast<std::size_t>(w)] == in_b)
                {
                    stamp_[static_cast<std::size_t>(w)] = seen;
                    const auto dw = static_cast<std::uint32_t>(out_adj_[w].size());
                    dfs_stack_.push_back({w, 0, dw ? static_cast<std::uint32_t>(rng.below(dw)) : 0});
                    pushed = true;
                    break;
                }
            }
            if (!pushed)
            {
                post_.push_back(dfs_stack_.back().v);
                dfs_stack_.pop_back();
            }
        }
    }
    // reverse post-order
    auto i = b;
    for (auto it = post_.rbegin(); it != post_.rend(); ++it, ++i)
    {
        order_[static_cast<std::size_t>(i)] = *it;
        pos_[static_cast<std::size_t>(*it)] = i;
    }
    touch(b);
}

void PseudoTopoOrder::order_free(std::int32_t b, std::int32_t c, Vertex anchor, const OpenerOptions& opt, Rng& rng,
                                 OpenStats& st)
{
    if (opt.border == BOrderRule::ShuffleOnly || c - b < 2)
        return;
    if (opt.border == BOrderRule::DfsPostorder)
    {
        order_free_dfs(b, c, anchor, rng);
        st.sort_ops += c - b;
        return;
    }
    if (opt.sweep)
        for (std::int32_t r = c - 1; r > b; --r)
        {
            const auto j = earlier_out_neighbour(b, r, opt.border, rng);
            ++st.sort_ops;
            if (j != -1)
                swap_positions(j, r);
        }
    std::int64_t iters = opt.sort_iters;
    if (opt.sort_iters_per_free > 0)
        iters = static_cast<std::int64_t>(std::ceil(opt.sort_iters_per_free * (c - b)));
    for (std::int64_t k = 0; k < iters; ++k)
    {
        const auto r = b + static_cast<std::int32_t>(rng.below(static_cast<std::uint64_t>(c - b)));
        const auto j = earlier_out_neighbour(b, r, opt.border, rng);
        ++st.sort_ops;
        if (j != -1)
            swap_positions(j, r);
    }
}

bool PseudoTopoOrder::open_pass(const OpenerOptions& opt, Rng& rng, const Deadline& deadline, OpenStats& st)
{
    ++st.passes;
    st.last_pass_sideways = false;
    const Weight old = value();
    std::vector<std::int32_t> pp = path_positions();
    // mark the current path's vertices (to detect plateau moves to a different vertex set)
    std::uint32_t path_mark = 0;
    if (opt.sideways > 0)
    {
        if (stamp_.size() != static_cast<std::size_t>(n_))
            stamp_.assign(static_cast<std::size_t>(n_), 0);
        if (cur_stamp_ > 0x7ffffff0U)
        {
            std::ranges::fill(stamp_, 0);
            cur_stamp_ = 0;
        }
        path_mark = 2 * (++cur_stamp_);
        for (auto i : pp)
            stamp_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])] = path_mark;
    }
    auto plateau_move = [&]() -> bool {
        if (opt.sideways <= 0 || value() != old || rng.uniform01() >= opt.sideways)
            return false;
        for (auto i : path_positions())
            if (stamp_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])] != path_mark)
            {
                ++st.sideways_moves;
                st.last_pass_sideways = true;
                return true;
            }
        return false;
    };
    if (opt.free_first && !pp.empty())
        pp.erase(pp.begin()); // 2018 FillPath never stores the first vertex of the DP path
    for (std::int32_t comp = 0; comp < scc_.num_components; ++comp)
    {
        const std::int32_t a = scc_.start[static_cast<std::size_t>(comp)];
        const std::int32_t d = scc_.start[static_cast<std::size_t>(comp) + 1];
        if (d - a <= opt.skip_block_size || d - a < 2)
            continue;
        // Move the path vertices of this block to its front (they are contiguous in the path,
        // and their positions are increasing).
        auto it = std::ranges::lower_bound(pp, a);
        // last path vertex before this block (anchor of the gap h = 0)
        const Vertex before_block = it == pp.begin() ? kNoVertex : order_[static_cast<std::size_t>(*(it - 1))];
        std::int32_t aa = a;
        for (; it != pp.end() && *it < d; ++it, ++aa)
            if (*it != aa)
                swap_positions(aa, *it); // aa < *it; later path positions are untouched
        std::int32_t b = aa;
        std::int32_t c = d;
        if (c - b < opt.min_free || c - b < 1)
            continue;
        const std::int32_t total = b - a; // number of path vertices in the block
        std::vector<std::int32_t> gaps(static_cast<std::size_t>(total) + 1);
        std::iota(gaps.begin(), gaps.end(), 0);
        if (opt.random_gap_order)
            rng.shuffle(gaps.begin(), gaps.end());
        // shuffle the free vertices once (2018)
        rng.shuffle(order_.begin() + b, order_.begin() + c);
        for (auto i = b; i < c; ++i)
            pos_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])] = i;
        touch(b);
        for (const std::int32_t h : gaps)
        {
            if (deadline.expired())
                return false;
            transfer(a, b, c, h);
            b = a + h;
            c = d - (total - h);
            if (opt.reshuffle_each_gap)
            {
                rng.shuffle(order_.begin() + b, order_.begin() + c);
                for (auto i = b; i < c; ++i)
                    pos_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])] = i;
                touch(b);
            }
            const Vertex anchor = h > 0 ? order_[static_cast<std::size_t>(a + h - 1)] : before_block;
            order_free(b, c, anchor, opt, rng, st);
            ++st.openings;
            evaluate(&st);
            if (value() > old)
            {
                ++st.improvements;
                return true;
            }
            if (plateau_move())
                return false;
            if (opt.try_reverse)
            {
                std::reverse(order_.begin() + b, order_.begin() + c);
                for (auto i = b; i < c; ++i)
                    pos_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])] = i;
                touch(b);
                ++st.openings;
                evaluate(&st);
                if (value() > old)
                {
                    ++st.improvements;
                    return true;
                }
            }
        }
    }
    return false;
}

} // namespace lsp
