#include "lsp/dfs.hpp"

#include <stdexcept>

namespace lsp {

InPlaceDfs::InPlaceDfs(Vertex n, const OrderedAdjacency& adj) : adj_(adj), used_(static_cast<std::size_t>(n), 0)
{
    seq_.reserve(static_cast<std::size_t>(n));
}

void InPlaceDfs::reset(const std::vector<Vertex>& seq)
{
    for (const auto& e : seq_)
        used_[static_cast<std::size_t>(e.v)] = 0;
    seq_.clear();
    value_ = 0;
    for (std::size_t i = 0; i < seq.size(); ++i)
    {
        const Vertex v = seq[i];
        if (used_[static_cast<std::size_t>(v)])
            throw std::logic_error("InPlaceDfs::reset: repeated vertex");
        if (i == 0)
        {
            seq_.push_back({v, -1, 0});
            used_[static_cast<std::size_t>(v)] = 1;
            continue;
        }
        const auto nb = adj_[seq[i - 1]];
        std::int32_t idx = -1;
        for (std::size_t j = 0; j < nb.size(); ++j)
            if (nb[j].v == v)
            {
                idx = static_cast<std::int32_t>(j);
                break;
            }
        if (idx < 0)
            throw std::logic_error("InPlaceDfs::reset: not a walk of the adjacency");
        push(seq[i - 1], idx);
    }
    best_ = seq_;
    best_value_ = value_;
    best_dirty_ = false;
    common_ = seq_.size();
}

std::int32_t InPlaceDfs::first_free(Vertex u, std::int32_t from) const
{
    const auto nb = adj_[u];
    for (auto i = static_cast<std::size_t>(from); i < nb.size(); ++i)
        if (!used_[static_cast<std::size_t>(nb[i].v)])
            return static_cast<std::int32_t>(i);
    return -1;
}

void InPlaceDfs::push(Vertex u, std::int32_t i)
{
    const Arc& a = adj_[u][static_cast<std::size_t>(i)];
    seq_.push_back({a.v, i, a.w});
    used_[static_cast<std::size_t>(a.v)] = 1;
    value_ += a.w;
}

void InPlaceDfs::pop()
{
    if (best_dirty_)
        sync_best();
    const Entry e = seq_.back();
    seq_.pop_back();
    used_[static_cast<std::size_t>(e.v)] = 0;
    value_ -= e.w;
    common_ = std::min(common_, seq_.size());
}

void InPlaceDfs::sync_best()
{
    best_.resize(common_);
    best_.insert(best_.end(), seq_.begin() + static_cast<std::ptrdiff_t>(common_), seq_.end());
    common_ = seq_.size();
    best_dirty_ = false;
}

void InPlaceDfs::greedy_extend()
{
    while (true)
    {
        const Vertex u = seq_.back().v;
        const std::int32_t i = first_free(u, 0);
        if (i < 0)
            return;
        push(u, i);
    }
}

bool InPlaceDfs::step()
{
    Vertex u = seq_.back().v;
    std::int32_t i = first_free(u, 0);
    while (i < 0 && seq_.size() > 1)
    {
        const std::int32_t popped_idx = seq_.back().idx;
        pop();
        u = seq_.back().v;
        i = first_free(u, popped_idx + 1);
    }
    if (i < 0)
        return false;
    push(u, i);
    greedy_extend();
    return true;
}

std::vector<Vertex> InPlaceDfs::sequence() const
{
    std::vector<Vertex> r;
    r.reserve(seq_.size());
    for (const auto& e : seq_)
        r.push_back(e.v);
    return r;
}

std::vector<Vertex> InPlaceDfs::run(const StopRule& stop, const Deadline& deadline, DfsStats& stats, const ValueSink& sink)
{
    if (value_ > best_value_)
    {
        best_value_ = value_;
        best_dirty_ = true;
        if (sink)
            sink(best_value_);
    }
    Stopwatch since;
    std::int64_t steps_since = 0;
    std::int64_t local = 0;
    const bool use_time = stop.seconds_without_improvement > 0;
    const bool use_steps = stop.steps_without_improvement > 0;
    while (true)
    {
        if ((local & 7) == 0)
        {
            if (deadline.expired())
                break;
            if (use_time && since.seconds() >= stop.seconds_without_improvement)
                break;
        }
        if (use_steps && steps_since >= stop.steps_without_improvement)
            break;
        if (!step())
            break;
        ++local;
        ++stats.steps;
        ++steps_since;
        if (value_ > best_value_)
        {
            best_value_ = value_;
            best_dirty_ = true;
            ++stats.improvements;
            steps_since = 0;
            if (use_time)
                since.reset();
            if (sink)
                sink(best_value_);
        }
    }
    if (best_dirty_)
        sync_best();
    std::vector<Vertex> r;
    r.reserve(best_.size());
    for (const auto& e : best_)
        r.push_back(e.v);
    return r;
}

std::vector<Vertex> reversed(std::vector<Vertex> v)
{
    std::ranges::reverse(v);
    return v;
}

} // namespace lsp
