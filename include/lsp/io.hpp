// Readers and writers for the graph formats used by the benchmarks.
#pragma once

#include "lsp/graph.hpp"

#include <iosfwd>
#include <string>

namespace lsp {

enum class GraphFormat
{
    Auto,     // decide from the file extension / content
    EdgeList, // "n m" header, then m lines "u v [w]", 0-based (format of code_2018)
    AdjList,  // line i (1-based) = out-neighbours of vertex i, 1-based (Farnebaeck's GoGameGraphs)
    Snap,     // lines "u v [w]" with '#' comments, arbitrary integer ids (SNAP); ids remapped, kept as labels
    Dimacs,   // "p sp n m", "a u v w", 1-based (9th DIMACS challenge)
};

GraphFormat parse_format(const std::string& name);
Digraph read_graph(const std::string& filename, GraphFormat format = GraphFormat::Auto);
Digraph read_graph(std::istream& in, GraphFormat format);

void write_edgelist(std::ostream& out, const Digraph& g);

// Path files: whitespace separated vertex indices (0-based, internal numbering).
std::vector<Vertex> read_path(const std::string& filename);
void write_path(std::ostream& out, const std::vector<Vertex>& p);

} // namespace lsp
