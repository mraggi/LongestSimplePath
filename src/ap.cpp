#include "lsp/ap.hpp"

#include <algorithm>
#include <functional>
#include <stdexcept>

namespace lsp {

namespace {

using I64 = std::int64_t;
constexpr I64 kMinus = std::numeric_limits<I64>::min() / 4;
__extension__ using I128 = __int128;

} // namespace

ApSolution solve_ap(const Digraph& g, Vertex source, const Deadline& deadline)
{
    ApSolution res;
    const Vertex n = g.num_vertices();
    if (n == 0)
    {
        res.solved = res.certified = true;
        return res;
    }
    const bool rooted = source != kNoVertex;
    if (rooted && (source < 0 || source >= n))
        throw std::invalid_argument("solve_ap: source out of range");
    const auto N = static_cast<std::size_t>(n);
    const std::int32_t D = n; // the dummy vertex, as a column (and as row n when there is no source)
    // Rows: the vertices, plus D without a source. With a source, D->s is forced, so the D row and
    // the column s are removed (no row lists s). Columns: the vertices and D.
    const std::int32_t rows = rooted ? n : n + 1;
    const I64 scale = static_cast<I64>(rows) + 1; // epsilon = 1 is then optimal (total slack < one unit)

    Weight wabs = 0;
    for (Vertex u = 0; u < n; ++u)
        for (const auto& a : g.out(u))
            wabs = std::max(wabs, a.w < 0 ? -a.w : a.w);
    if (wabs > (I64{1} << 40) / scale)
        throw std::overflow_error("solve_ap: arc weights too large for the scaled auction");
    const I64 cmax = std::max<I64>(wabs * scale, 1);

    // rows 0..n-1 in CSR form: (column, scaled benefit). The D row (benefit 0 on every column) is implicit.
    std::vector<std::size_t> off(N + 1, 0);
    std::vector<std::int32_t> col;
    std::vector<I64> ben;
    col.reserve(g.num_edges() + 2 * N);
    ben.reserve(g.num_edges() + 2 * N);
    for (Vertex u = 0; u < n; ++u)
    {
        if (!rooted || u != source)
        {
            col.push_back(u); // u skipped
            ben.push_back(0);
        }
        col.push_back(D); // u ends the path
        ben.push_back(0);
        for (const auto& a : g.out(u))
        {
            if (rooted && a.v == source)
                continue;
            col.push_back(a.v);
            ben.push_back(a.w * scale);
        }
        off[static_cast<std::size_t>(u) + 1] = col.size();
    }

    std::vector<I64> price(N + 1, 0);
    std::vector<std::int32_t> owner(N + 1, -1);
    std::vector<std::int32_t> assigned(static_cast<std::size_t>(rows), -1);

    // Lazy min-heap of (price, column) for the D row. Prices only increase, so an entry is current
    // iff its price equals the column's price.
    using HeapEntry = std::pair<I64, std::int32_t>;
    std::vector<HeapEntry> heap;
    const auto cmp = std::greater<HeapEntry>{};
    auto heap_rebuild = [&] {
        heap.clear();
        for (std::int32_t j = 0; j <= n; ++j)
            heap.emplace_back(price[static_cast<std::size_t>(j)], j);
        std::make_heap(heap.begin(), heap.end(), cmp);
    };
    auto heap_push = [&](std::int32_t j) {
        heap.emplace_back(price[static_cast<std::size_t>(j)], j);
        std::push_heap(heap.begin(), heap.end(), cmp);
        if (heap.size() > 4 * N + 64)
            heap_rebuild();
    };
    auto heap_drop_stale = [&] {
        while (!heap.empty() && heap.front().first != price[static_cast<std::size_t>(heap.front().second)])
        {
            std::pop_heap(heap.begin(), heap.end(), cmp);
            heap.pop_back();
        }
    };
    if (!rooted)
        heap_rebuild();

    // epsilon-scaling forward auction (Gauss-Seidel bids, FIFO queue of unassigned rows)
    const auto R = static_cast<std::size_t>(rows);
    std::vector<std::int32_t> queue(R);
    constexpr I64 theta = 5;
    I64 eps = std::max<I64>(1, cmax / theta);
    for (;;)
    {
        ++res.phases;
        std::ranges::fill(owner, -1);
        std::ranges::fill(assigned, -1);
        for (std::size_t i = 0; i < R; ++i)
            queue[i] = static_cast<std::int32_t>(i);
        std::size_t qh = 0;
        std::size_t qn = R;
        while (qn > 0)
        {
            const std::int32_t i = queue[qh];
            qh = qh + 1 == R ? 0 : qh + 1;
            --qn;
            if ((++res.bids & 4095) == 0 && deadline.expired())
                return res;
            std::int32_t j1 = -1;
            I64 v1 = kMinus;
            I64 v2 = kMinus;
            if (i == n) // the D row: benefit 0 on every column, so the best column has the lowest price
            {
                heap_drop_stale();
                j1 = heap.front().second;
                v1 = -heap.front().first;
                std::pop_heap(heap.begin(), heap.end(), cmp);
                heap.pop_back();
                heap_drop_stale();
                if (!heap.empty())
                    v2 = -heap.front().first;
            }
            else
            {
                const auto iu = static_cast<std::size_t>(i);
                for (std::size_t k = off[iu]; k < off[iu + 1]; ++k)
                {
                    const I64 v = ben[k] - price[static_cast<std::size_t>(col[k])];
                    if (v > v1)
                    {
                        v2 = v1;
                        v1 = v;
                        j1 = col[k];
                    }
                    else if (v > v2)
                        v2 = v;
                }
            }
            if (v2 == kMinus) // a single admissible column: any increment keeps epsilon-CS
                v2 = v1 - cmax - eps;
            const auto ju = static_cast<std::size_t>(j1);
            price[ju] += (v1 - v2) + eps;
            if (!rooted)
                heap_push(j1);
            if (owner[ju] >= 0)
            {
                assigned[static_cast<std::size_t>(owner[ju])] = -1;
                queue[(qh + qn) % R] = owner[ju];
                ++qn;
            }
            owner[ju] = i;
            assigned[static_cast<std::size_t>(i)] = j1;
        }
        if (eps == 1)
            break;
        eps = std::max<I64>(1, eps / theta);
    }
    res.solved = true;

    // value (unscaled, from the digraph) and the dual certificate
    for (Vertex u = 0; u < n; ++u)
    {
        const std::int32_t j = assigned[static_cast<std::size_t>(u)];
        if (j != D && j != u)
            res.bound += *g.arc_weight(u, j);
    }
    I128 dual = 0;
    for (Vertex u = 0; u < n; ++u)
    {
        I64 best = kMinus;
        const auto iu = static_cast<std::size_t>(u);
        for (std::size_t k = off[iu]; k < off[iu + 1]; ++k)
            best = std::max(best, ben[k] - price[static_cast<std::size_t>(col[k])]);
        dual += best;
    }
    if (!rooted)
        dual -= *std::ranges::min_element(price); // the D row
    for (std::int32_t j = 0; j <= n; ++j)
        if (!(rooted && j == source))
            dual += price[static_cast<std::size_t>(j)];
    res.certified = dual - static_cast<I128>(res.bound) * scale < scale;

    // the D-cycle and the other cycles
    std::vector<char> seen(N, 0);
    Vertex v = rooted ? source : assigned[N];
    while (v != D)
    {
        res.path.push_back(v);
        seen[static_cast<std::size_t>(v)] = 1;
        v = assigned[static_cast<std::size_t>(v)];
    }
    for (Vertex s = 0; s < n; ++s)
    {
        if (seen[static_cast<std::size_t>(s)] || assigned[static_cast<std::size_t>(s)] == s)
            continue;
        std::vector<Vertex> c;
        for (v = s; !seen[static_cast<std::size_t>(v)]; v = assigned[static_cast<std::size_t>(v)])
        {
            seen[static_cast<std::size_t>(v)] = 1;
            c.push_back(v);
        }
        res.cycles.push_back(std::move(c));
    }
    return res;
}

std::vector<Vertex> patch_ap(const Digraph& g, const ApSolution& ap, Vertex source, PatchStats* stats)
{
    const Vertex n = g.num_vertices();
    const auto N = static_cast<std::size_t>(n);
    const bool rooted = source != kNoVertex;
    PatchStats st;
    auto W = [&](Vertex a, Vertex b) { return *g.arc_weight(a, b); };
    auto at = [](auto& vec, Vertex v) -> auto& { return vec[static_cast<std::size_t>(v)]; };

    // the path as a doubly linked list; nw[v] = weight of the arc v -> nxt[v]
    std::vector<Vertex> nxt(N, kNoVertex);
    std::vector<Vertex> prv(N, kNoVertex);
    std::vector<Weight> nw(N, 0);
    std::vector<char> inP(N, 0);
    Vertex head = kNoVertex;
    Vertex tail = kNoVertex;
    for (Vertex v : ap.path)
    {
        if (tail == kNoVertex)
            head = v;
        else
        {
            at(nxt, tail) = v;
            at(prv, v) = tail;
            at(nw, tail) = W(tail, v);
        }
        tail = v;
        at(inP, v) = 1;
    }

    // cycles: cin[k][r] = weight of the cycle arc entering C[r]
    const auto& cycles = ap.cycles;
    const auto K = static_cast<std::int32_t>(cycles.size());
    std::vector<std::int32_t> cyc(N, -1);
    std::vector<std::int32_t> cidx(N, -1);
    std::vector<std::vector<Weight>> cin(cycles.size());
    std::vector<Weight> cw(cycles.size(), 0);
    std::vector<char> alive(cycles.size(), 1);
    for (std::int32_t k = 0; k < K; ++k)
    {
        const auto& C = cycles[static_cast<std::size_t>(k)];
        const auto L = C.size();
        auto& ci = cin[static_cast<std::size_t>(k)];
        ci.resize(L);
        for (std::size_t r = 0; r < L; ++r)
        {
            at(cyc, C[r]) = k;
            at(cidx, C[r]) = static_cast<std::int32_t>(r);
            ci[r] = W(C[(r + L - 1) % L], C[r]);
            cw[static_cast<std::size_t>(k)] += ci[r];
        }
    }

    // Puts the rotation of cycle k that starts at C[r] (and ends at C[r-1]) between a and b;
    // a = kNoVertex: at the start; b = kNoVertex: at the end.
    auto splice = [&](Vertex a, Vertex b, std::int32_t k, std::size_t r) {
        const auto& C = cycles[static_cast<std::size_t>(k)];
        const auto L = C.size();
        Vertex prev = a;
        for (std::size_t t = 0; t < L; ++t)
        {
            const Vertex v = C[(r + t) % L];
            at(inP, v) = 1;
            at(cyc, v) = -1;
            if (prev == kNoVertex)
            {
                head = v;
                at(prv, v) = kNoVertex;
            }
            else
            {
                at(nxt, prev) = v;
                at(prv, v) = prev;
                at(nw, prev) = t == 0 ? W(prev, v) : cin[static_cast<std::size_t>(k)][(r + t) % L];
            }
            prev = v;
        }
        if (b == kNoVertex)
        {
            tail = prev;
            at(nxt, prev) = kNoVertex;
        }
        else
        {
            at(nxt, prev) = b;
            at(prv, b) = prev;
            at(nw, prev) = W(prev, b);
        }
        alive[static_cast<std::size_t>(k)] = 0;
    };

    enum class Move { None, Append, Prepend, Insert, CutSuffix, CutPrefix };

    // (1) merges with a gain
    auto gainful = [&] {
        bool changed = true;
        while (changed)
        {
            changed = false;
            for (std::int32_t k = 0; k < K; ++k)
            {
                const auto ku = static_cast<std::size_t>(k);
                if (!alive[ku])
                    continue;
                const auto& C = cycles[ku];
                const auto& ci = cin[ku];
                const auto L = C.size();
                if (head == kNoVertex) // empty path (no source): open the cycle at its lightest arc
                {
                    const auto r = static_cast<std::size_t>(std::ranges::min_element(ci) - ci.begin());
                    splice(kNoVertex, kNoVertex, k, r);
                    ++st.merged;
                    changed = true;
                    continue;
                }
                Weight best = 0;
                Move move = Move::None;
                std::size_t br = 0;
                Vertex ba = kNoVertex;
                for (const auto& a : g.out(tail)) // append: tail -> C[r] ... C[r-1]
                    if (at(cyc, a.v) == k)
                    {
                        const auto r = static_cast<std::size_t>(at(cidx, a.v));
                        const Weight gain = a.w + cw[ku] - ci[r];
                        if (gain > best)
                            best = gain, move = Move::Append, br = r;
                    }
                if (!rooted)
                    for (const auto& a : g.in(head)) // prepend: C[r] ... C[r-1] = a.v -> head
                        if (at(cyc, a.v) == k)
                        {
                            const auto r = (static_cast<std::size_t>(at(cidx, a.v)) + 1) % L;
                            const Weight gain = cw[ku] - ci[r] + a.w;
                            if (gain > best)
                                best = gain, move = Move::Prepend, br = r;
                        }
                for (std::size_t r = 0; r < L; ++r) // insert between a and a' = nxt[a]: a -> C[r] ... C[r-1] -> a'
                {
                    const Vertex c = C[(r + L - 1) % L];
                    const Weight seg = cw[ku] - ci[r];
                    for (const auto& a : g.in(C[r]))
                    {
                        if (!at(inP, a.v) || a.v == tail)
                            continue;
                        const auto wb = g.arc_weight(c, at(nxt, a.v));
                        if (!wb)
                            continue;
                        const Weight gain = a.w + seg + *wb - at(nw, a.v);
                        if (gain > best)
                            best = gain, move = Move::Insert, br = r, ba = a.v;
                    }
                }
                if (move == Move::None)
                    continue;
                if (move == Move::Append)
                    splice(tail, kNoVertex, k, br);
                else if (move == Move::Prepend)
                    splice(kNoVertex, head, k, br);
                else
                    splice(ba, at(nxt, ba), k, br);
                ++st.merged;
                changed = true;
            }
        }
    };

    // (2) the best merge that cuts off a suffix or a prefix of the path; false if none gains
    std::vector<Weight> pref(N, 0); // pref[v] = weight of the path from head to v
    auto lossy = [&]() -> bool {
        if (head == kNoVertex)
            return false;
        Weight tot = 0;
        for (Vertex v = head; v != kNoVertex; v = at(nxt, v))
        {
            at(pref, v) = tot;
            if (at(nxt, v) != kNoVertex)
                tot += at(nw, v);
        }
        Weight best = 0;
        Move move = Move::None;
        std::int32_t bk = -1;
        std::size_t br = 0;
        Vertex bv = kNoVertex;
        for (std::int32_t k = 0; k < K; ++k)
        {
            const auto ku = static_cast<std::size_t>(k);
            if (!alive[ku])
                continue;
            const auto& C = cycles[ku];
            const auto L = C.size();
            for (std::size_t r = 0; r < L; ++r)
            {
                const Weight seg = cw[ku] - cin[ku][r];
                for (const auto& a : g.in(C[r])) // path up to a.v, then a.v -> C[r] ... C[r-1]
                    if (at(inP, a.v))
                    {
                        const Weight gain = at(pref, a.v) + a.w + seg - tot;
                        if (gain > best)
                            best = gain, move = Move::CutSuffix, bk = k, br = r, bv = a.v;
                    }
                if (!rooted)
                    for (const auto& a : g.out(C[(r + L - 1) % L])) // C[r] ... C[r-1] -> a.v, then the path from a.v
                        if (at(inP, a.v))
                        {
                            const Weight gain = seg + a.w - at(pref, a.v);
                            if (gain > best)
                                best = gain, move = Move::CutPrefix, bk = k, br = r, bv = a.v;
                        }
            }
        }
        if (move == Move::None)
            return false;
        auto release = [&](Vertex v) {
            at(inP, v) = 0;
            at(nxt, v) = at(prv, v) = kNoVertex;
        };
        if (move == Move::CutSuffix)
        {
            for (Vertex v = at(nxt, bv); v != kNoVertex;)
            {
                const Vertex x = at(nxt, v);
                release(v);
                v = x;
            }
            at(nxt, bv) = kNoVertex;
            tail = bv;
            splice(bv, kNoVertex, bk, br);
        }
        else
        {
            for (Vertex v = at(prv, bv); v != kNoVertex;)
            {
                const Vertex x = at(prv, v);
                release(v);
                v = x;
            }
            at(prv, bv) = kNoVertex;
            head = bv;
            splice(kNoVertex, bv, bk, br);
        }
        ++st.lossy;
        return true;
    };

    gainful();
    while (lossy())
        gainful();

    // (3) greedy extension of both ends
    if (head == kNoVertex)
    {
        if (!rooted)
        {
            if (stats)
                *stats = st;
            return {};
        }
        head = tail = source;
        at(inP, source) = 1;
    }
    for (bool ext = true; ext;)
    {
        ext = false;
        Vertex best = kNoVertex;
        Weight bw = 0;
        std::size_t bd = 0;
        for (const auto& a : g.out(tail))
            if (!at(inP, a.v) && a.w > 0 && (best == kNoVertex || a.w > bw || (a.w == bw && g.outdegree(a.v) > bd)))
                best = a.v, bw = a.w, bd = g.outdegree(a.v);
        if (best != kNoVertex)
        {
            at(nxt, tail) = best;
            at(prv, best) = tail;
            at(nw, tail) = bw;
            at(nxt, best) = kNoVertex;
            at(inP, best) = 1;
            tail = best;
            ++st.extended;
            ext = true;
        }
        if (rooted)
            continue;
        best = kNoVertex;
        for (const auto& a : g.in(head))
            if (!at(inP, a.v) && a.w > 0 && (best == kNoVertex || a.w > bw || (a.w == bw && g.indegree(a.v) > bd)))
                best = a.v, bw = a.w, bd = g.indegree(a.v);
        if (best != kNoVertex)
        {
            at(prv, head) = best;
            at(nxt, best) = head;
            at(nw, best) = bw;
            at(prv, best) = kNoVertex;
            at(inP, best) = 1;
            head = best;
            ++st.extended;
            ext = true;
        }
    }

    std::vector<Vertex> p;
    for (Vertex v = head; v != kNoVertex; v = at(nxt, v))
        p.push_back(v);
    for (std::int32_t k = 0; k < K; ++k)
        st.left += alive[static_cast<std::size_t>(k)] ? 1 : 0;
    if (stats)
        *stats = st;
    return p;
}

} // namespace lsp
