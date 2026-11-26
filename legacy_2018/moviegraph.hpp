#pragma once
#include "digraph.hpp"
class Path;

class MovieGraph : public DiGraph
{
public:
    MovieGraph(const std::vector< std::string >& mvnames);
    void PrintPath(const Path& P);
private:
    // Utils for creating the graph
	void add_links(size_t i);
	std::vector<node_t> movies_that_start_with(const std::string& name);
// 	void branch_and_bound();
	
};

bool SecondStartsWithFirst(const std::string& first, const std::string& second);
