#include "genomegraph.hpp"
#include "path.hpp"
#include "pseudotopoorder.hpp"
#include <sstream>
#include <unordered_set>

template <class T>
std::vector<T> sorted(std::vector<T> V)
{
    sort(V.begin(), V.end());
    return V;
}

GenomeGraph::GenomeGraph(const std::vector< std::string >& mvnames) : 	DiGraph(mvnames)
{
	for (node_t i = 0; i < m_n; ++i)
	{
// 		std::cout << "vertex values " << i << std::endl;
// 		std::cout << "Adding links " << i << std::endl;
		add_links(i);
// 		std::cout << "Done!" << std::endl;
	}
}

void GenomeGraph::add_links(size_t i)
{
	std::unordered_set<node_t> alreadyadded;
	const std::string& currgenome = get_vertex_name(i);
	int length = currgenome.size();
// 	std::cout << "Adding links to " << currgenome << " # = " << length <<std::endl;
	std::cout << "Completed: " << double(i*100)/num_vertices() << "%" << std::endl;
	for (int l = MIN_SIZE_INTERSECTION; l < length; ++l)
	{
		std::string cgenomeend(currgenome.begin()+(length-l),currgenome.end());
		auto V = genomes_that_start_with(cgenomeend);
		for (auto v : V)
		{
			if (v != i && alreadyadded.count(v) == 0)
			{
				weight_t weight = l;
				add_edge(i,v,weight);
				alreadyadded.insert(v);
// 				std::cout << "\tadding edge: " << i << " --> " << v << " #= " << weight << std::endl;
			}
		}
	}
	
}


std::vector<node_t> GenomeGraph::genomes_that_start_with(const std::string& name)
{
// 	std::cout << "name = " << name << std::endl;
	std::vector<node_t> toReturn;
	
	// true true true false false (regresa el primer false)
	auto it = partition_point(get_vertex_names().begin(), get_vertex_names().end(), [&name](const std::string& x)
	{
		return x < name;
	});
	
	while(it != get_vertex_names().end() && SecondStartsWithFirst(name,*it) )
	{
		if (name.size() < it->size())
		{
			toReturn.push_back(it-get_vertex_names().begin());
		}
		++it;
	}
	return toReturn;
}

// bool SecondStartsWithFirst(const std::string& first, const std::string& second)
// {
// 	int fs = first.size();
// 	int ss = second.size();
// 	if (fs > ss)
// 		return false;
// 	
// 	for (size_t i = 0; i < fs; ++i)
// 	{
// 		if (first[i] != second[i])
// 			return false;
// 	}
// 	
// 	return true;
// }
/*
void GenomeGraph::PrintPath(const Path& P)
{
    auto A = P.get_path();
//     reverse(A.begin(), A.end());
	std::cout << get_vertex_name(A.front());
	auto i = A.begin();
	node_t prevnode = *i;
	++i;
    for ( ; i != A.end(); ++i)
    {
		node_t currnode = *i;
        edge_value(prevnode,currnode);
		prevnode = currnode;
    }
    std::cout << std::endl;
}*/
