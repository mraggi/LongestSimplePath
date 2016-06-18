#include "moviegraph.hpp"
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

MovieGraph::MovieGraph(const vector< string >& mvnames) : 	DiGraph(mvnames)
{
	for (node_t i = 0; i < m_n; ++i)
	{
// 		cout << "vertex values " << i << endl;
// 		cout << "Adding links " << i << endl;
		add_links(i);
// 		cout << "Done!" << endl;
	}
}

void MovieGraph::add_links(size_t i)
{
	unordered_set<node_t> alreadyadded;
	const string& movie = get_vertex_name(i);
	size_t length = movie.size();
// 	cout << "Adding links to " << movie << " # = " << length <<endl;
	for (int l = length-1; l > 0; --l)
	{
		if (movie[l] != ' ')
		{
// 			spacefound = true;
			continue;
		}
// 		if (!spacefound)
// 			continue;
		string movieend(movie.begin()+(l+1),movie.end());
		auto V =  movies_that_start_with(movieend);
		for (auto v : V)
		{
			if (v != i && alreadyadded.count(v) == 0)
			{
				weight_t weight = (length+get_vertex_name(v).size())-movieend.size()*2;
				add_edge(i,v,weight);
				alreadyadded.insert(v);
// 				cout << "\tadding edge: " << movienames[v] << " #= " << float(movienames[v].size()) << " -> with weight " << float(weight) << endl;
			}
		}
	}
	
}


vector<node_t> MovieGraph::movies_that_start_with(const string& name)
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

bool SecondStartsWithFirst(const string& first, const string& second)
{
	int fs = first.size();
	int ss = second.size();
	if (fs > ss)
		return false;
	
	for (size_t i = 0; i < fs; ++i)
	{
		if (first[i] != second[i])
			return false;
	}
	
	return true;
}

int intersection_size(const string& a, const string& b)
{
	size_t length = a.size();
// 	cout << "Adding links to " << movie << " # = " << length <<endl;
	for (int l = length-1; l > 0; --l)
	{
		if (a[l] != ' ')
		{
// 			spacefound = true;
			continue;
		}
// 		if (!spacefound)
// 			continue;
		string aend(a.begin()+(l+1),a.end());
		string bbegin(b.begin(),b.begin()+aend.size());
        if (aend == bbegin)
            return aend.size();
        
	}
	return 0;
}

void MovieGraph::PrintPath(const Path& P)
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
		const string& movie = get_vertex_name(currnode);
        
		cout << movie.substr(intersection_size(get_vertex_name(prevnode),get_vertex_name(currnode)));
//         cout << movie;
		prevnode = currnode;
    }
    cout << endl;
}
