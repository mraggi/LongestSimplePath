#include "genomegraph.hpp"
#include "path.hpp"
#include "pseudotopoorder.hpp"
#include <sstream>
#include <unordered_set>

template <class T>
vector<T> sorted(vector<T> V)
{
    sort(V.begin(), V.end());
    return V;
}

GenomeGraph::GenomeGraph(const vector< string >& mvnames) : 	DiGraph(mvnames)
{
	for (node_t i = 0; i < m_n; ++i)
	{
// 		cout << "vertex values " << i << endl;
// 		cout << "Adding links " << i << endl;
		add_links(i);
// 		cout << "Done!" << endl;
	}
}

void GenomeGraph::add_links(size_t i)
{
	unordered_set<node_t> alreadyadded;
	const string& currgenome = get_vertex_name(i);
	int length = currgenome.size();
// 	cout << "Adding links to " << currgenome << " # = " << length <<endl;
	cout << "Completed: " << double(i*100)/num_vertices() << "%" << endl;
	for (int l = MIN_SIZE_INTERSECTION; l < length; ++l)
	{
		string cgenomeend(currgenome.begin()+(length-l),currgenome.end());
		auto V = genomes_that_start_with(cgenomeend);
		for (auto v : V)
		{
			if (v != i && alreadyadded.count(v) == 0)
			{
				weight_t weight = l;
				add_edge(i,v,weight);
				alreadyadded.insert(v);
// 				cout << "\tadding edge: " << i << " --> " << v << " #= " << weight << endl;
			}
		}
	}
	
}


vector<node_t> GenomeGraph::genomes_that_start_with(const string& name)
{
// 	cout << "name = " << name << endl;
	vector<node_t> toReturn;
	
	// true true true false false (regresa el primer false)
	auto it = partition_point(get_vertex_names().begin(), get_vertex_names().end(), [&name](const string& x)
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

// bool SecondStartsWithFirst(const string& first, const string& second)
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
	cout << get_vertex_name(A.front());
	auto i = A.begin();
	node_t prevnode = *i;
	++i;
    for ( ; i != A.end(); ++i)
    {
		node_t currnode = *i;
        edge_value(prevnode,currnode);
		prevnode = currnode;
    }
    cout << endl;
}*/
