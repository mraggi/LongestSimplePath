// The complete heuristic (2018 schedule and variants), baselines and the option set used for ablations.
#pragma once

#include "lsp/ap.hpp"
#include "lsp/dfs.hpp"
#include "lsp/pto.hpp"

#include <map>

namespace lsp {

enum class PtoMode
{
    Open,     // edge opening (paper §3.2.4-3.2.5)
    Reimpose, // only impose the path on fresh random PTOs (paper §3.2.1)
};

struct SolverOptions
{
    // A default-constructed SolverOptions is the recommended method (preset "newap") with a time
    // limit of kDefaultTimeLimit seconds. The member initializers below are the values of preset
    // "2018" (faithful to the 2016-2018 code); every preset starts from them (apply_preset).
    SolverOptions();
    static constexpr double kDefaultTimeLimit = 10;

    std::string method = "classic"; // classic | ils | rdfs | greedy | gamm
    double time_limit = 0;           // seconds (0 = no limit); kDefaultTimeLimit when default-constructed
    std::uint64_t seed = 1;
    Vertex source = kNoVertex;       // rooted variant: paths must start here
    bool largest_wcc_only = true;    // 2018 keeps only the largest weakly connected component

    // vertex orders / DFS guidance
    OrderingOptions ordering;
    int param_restarts = 10;                 // number of parameter sets (outer loop of 2018)
    std::string params = "2018";             // 2018 | fixed (always P0) | random
    // DFS
    int fwd_starts = 5;
    int bwd_starts = 3;
    int erase = 16;                          // vertices erased before searching in the other direction
    bool backward = true;                    // run dfs_backward_full
    bool other_direction = true;             // erase-and-extend in the other direction
    bool ranked_starts = true;               // alternate ranked and random starts (else random only)
    bool faithful_backward_greedy = true;    // 2018 extends the back greedily before a backward search
    StopRule dfs_stop{0.1, 0};
    // PTO
    bool pto = true;
    bool pto_each_dfs = true;
    PtoMode pto_mode = PtoMode::Open;
    StopRule pto_stop{4.0, 0};               // seconds / passes without improvement
    bool pto_incremental = true;
    OpenerOptions opener;
    std::int64_t sideways_max = 1000;        // plateau moves without improvement before a PTO search stops

    // ILS (modern schedule): see solver.cpp
    double ils_keep_fraction = 0.5;          // fraction of the path kept when perturbing
    std::int64_t ils_dfs_steps = 0;          // 0: use dfs_stop
    std::int64_t ils_rounds = 0;             // stop after this many ILS rounds (0: only the time limit); for deterministic runs
    std::string init_path;                   // ils: file with a starting path (original vertex ids); used if better than the DFS start
    bool ap_init = false;                    // ils: solve the assignment relaxation, patch its cycles into a path, start from it
                                             // if better than the DFS start; stop early if a path reaches the certified bound

    // baselines
    std::int64_t rdfs_steps = 0;             // restart after this many steps w/o improvement (0: 10*n)
    int ga_population = 100;                 // GAMM (Portugal et al. 2010)
    double ga_elite = 0.05;
    int ga_stagnation = 15;                  // generations without improvement before restarting

    // Presets: "newap" (the default: "new" plus the assignment-relaxation opening, ap_init=1),
    // "new" (anytime ILS, DFS post-order ordering of B, symmetric orders, work-based stopping
    // rules, no restriction to the largest weak component) and "2018" (faithful to the 2016-2018
    // code). A preset resets the algorithmic options and keeps time_limit, seed and source.
    void apply_preset(const std::string& name);
    // Applies "key=value"; throws on unknown keys.
    void set(const std::string& key, const std::string& value);
    [[nodiscard]] std::map<std::string, std::string> describe() const;

private:
    struct Base2018
    {
    };
    explicit SolverOptions(Base2018) {} // the member initializers: preset "2018"
};

struct TrajectoryPoint
{
    double time;
    Weight value;
};

// The assignment-relaxation opening (ap_init); see ap.hpp.
struct ApReport
{
    bool ran = false;
    bool solved = false;            // false if the deadline expired during the auction
    bool certified = false;         // dual check of the AP optimum passed
    bool bound_global = false;      // the searched graph is the whole problem (no largest-WCC restriction)
    Weight bound = 0;               // AP optimum: upper bound on every simple path of the searched graph
    std::size_t path_vertices = 0;  // vertices on the AP path
    std::size_t cycles = 0;         // AP cycles of length >= 2
    std::size_t cycle_vertices = 0;
    Weight patched = 0;             // value of the patched path
    double ap_seconds = 0;
    double patch_seconds = 0;
    std::int64_t bids = 0;
    int phases = 0;
    PatchStats patch;
};

struct SolveResult
{
    Path best;                      // in the numbering of the input digraph
    std::vector<TrajectoryPoint> trajectory; // (time since start of solve, value) at every improvement
    double preprocess_seconds = 0;
    double total_seconds = 0;
    DfsStats dfs;
    OpenStats open;
    std::int64_t pto_searches = 0;
    std::int64_t restarts = 0;
    Vertex work_vertices = 0;       // vertices of the graph actually searched
    std::size_t work_edges = 0;
    ApReport ap;
    bool proven_optimal = false;    // the best path reached the certified AP bound (of the whole problem)
};

SolveResult solve(const Digraph& g, const SolverOptions& opt);

} // namespace lsp
