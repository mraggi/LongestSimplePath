#pragma once

#include "digraph.hpp"

class Path;

const int MIN_SIZE_INTERSECTION = 6;

class GenomeGraph : public DiGraph
{
public:
    GenomeGraph(const std::vector<std::string>& genomes);
    void PrintPath(const Path& P);
private:
    // Utils for creating the graph
	void add_links(size_t i);
	std::vector<node_t> genomes_that_start_with(const std::string& name);
	
};

bool SecondStartsWithFirst(const std::string& first, const std::string& second);
