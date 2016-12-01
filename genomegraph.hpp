#pragma once

#include "digraph.hpp"

class Path;

const int MIN_SIZE_INTERSECTION = 6;

using namespace std;

class GenomeGraph : public DiGraph
{
public:
    GenomeGraph(const vector<string>& genomes);
    void PrintPath(const Path& P);
private:
    // Utils for creating the graph
	void add_links(size_t i);
	vector<node_t> genomes_that_start_with(const string& name);
	
};

bool SecondStartsWithFirst(const string& first, const string& second);
