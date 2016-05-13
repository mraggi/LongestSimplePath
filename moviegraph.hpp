#pragma once

#include "digraph.hpp"

class Path;

using namespace std;

class MovieGraph : public DiGraph
{
public:
    MovieGraph(const vector< string >& mvnames);
    void PrintPath(const Path& P);
private:
    // Utils for creating the graph
	void add_links(size_t i);
	vector<node_t> movies_that_start_with(const string& name);
// 	void branch_and_bound();
	
};

bool SecondStartsWithFirst(const string& first, const string& second);
