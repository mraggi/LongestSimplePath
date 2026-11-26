// Core types, deterministic RNG and timing utilities.
#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace lsp {

using Vertex = std::int32_t;
using Weight = std::int64_t;

inline constexpr Vertex kNoVertex = -1;
// "minus infinity" for DP values; far from overflow when weights are added.
inline constexpr Weight kNegInf = std::numeric_limits<Weight>::min() / 4;

// xoshiro256** seeded through splitmix64. Deterministic on every platform,
// unlike std::uniform_int_distribution & co.
class Rng
{
public:
    explicit Rng(std::uint64_t seed = 0) { reseed(seed); }

    void reseed(std::uint64_t seed)
    {
        for (auto& x : s_)
        {
            seed += 0x9e3779b97f4a7c15ULL;
            std::uint64_t z = seed;
            z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
            x = z ^ (z >> 31U);
        }
    }

    std::uint64_t next()
    {
        const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
        const std::uint64_t t = s_[1] << 17U;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return result;
    }

    // Uniform integer in [0, n), n > 0 (Lemire's nearly divisionless method).
    std::uint64_t below(std::uint64_t n)
    {
        __uint128_t m = static_cast<__uint128_t>(next()) * n;
        auto low = static_cast<std::uint64_t>(m);
        if (low < n)
        {
            const std::uint64_t threshold = (0 - n) % n;
            while (low < threshold)
            {
                m = static_cast<__uint128_t>(next()) * n;
                low = static_cast<std::uint64_t>(m);
            }
        }
        return static_cast<std::uint64_t>(m >> 64U);
    }

    // Uniform integer in [lo, hi] (inclusive).
    std::int64_t between(std::int64_t lo, std::int64_t hi)
    {
        return lo + static_cast<std::int64_t>(below(static_cast<std::uint64_t>(hi - lo + 1)));
    }

    double uniform01() { return static_cast<double>(next() >> 11U) * 0x1.0p-53; }

    template <class It>
    void shuffle(It first, It last)
    {
        const auto n = static_cast<std::uint64_t>(last - first);
        for (std::uint64_t i = n; i > 1; --i)
        {
            const std::uint64_t j = below(i);
            std::iter_swap(first + static_cast<std::ptrdiff_t>(i - 1), first + static_cast<std::ptrdiff_t>(j));
        }
    }

private:
    static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    std::uint64_t s_[4]{};
};

using Clock = std::chrono::steady_clock;

class Stopwatch
{
public:
    Stopwatch() : start_(Clock::now()) {}
    void reset() { start_ = Clock::now(); }
    [[nodiscard]] double seconds() const
    {
        return std::chrono::duration<double>(Clock::now() - start_).count();
    }

private:
    Clock::time_point start_;
};

// Absolute deadline; "never" is represented by time_point::max().
struct Deadline
{
    Clock::time_point when = Clock::time_point::max();

    static Deadline never() { return {}; }
    static Deadline in(double seconds)
    {
        if (seconds <= 0 || seconds > 1e9)
            return never();
        return {Clock::now() + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(seconds))};
    }
    [[nodiscard]] bool is_never() const { return when == Clock::time_point::max(); }
    [[nodiscard]] bool expired() const { return !is_never() && Clock::now() >= when; }
};

} // namespace lsp
