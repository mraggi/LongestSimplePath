#include "digraph.hpp"
#include "path.hpp"
#include "pseudotopoorder.hpp"
#include <sstream>

//Heuristics

DiGraph::DiGraph(node_t numNodes) :
			m_n(numNodes),
			m_outgraph(numNodes),
			m_ingraph(numNodes),
// 			m_edge_values(numNodes),
			m_processed(false),
			m_strongly_connected_components(),
			m_scc_coloring(),
			m_scc_rank_out(),
			m_scc_rank_in(),
			m_weak_coloring(),
			m_weakly_connected_components(),
			m_heuristic_out(),
			m_heuristic_in(),
			m_basic_topological_ordering(),
			m_basic_topological_ordering_in(),
			m_basic_topological_ordering_inverse(),
			m_basic_topological_ordering_inverse_in(),
			m_scc_big_components(),
			m_node_names(numNodes),
			m_namemap()
{
	for (node_t i = 0; i < numNodes; ++i)
    {
        m_node_names[i] = to_string(i);
        m_namemap[to_string(i)] = i;
    }
    
}

DiGraph::DiGraph(const vector<string>& vnames) :
			m_n(vnames.size()),
			m_outgraph(vnames.size()),
			m_ingraph(vnames.size()),
// 			m_edge_values(vnames.size()),
			m_processed(false),
			m_strongly_connected_components(),
			m_scc_coloring(),
			m_scc_rank_out(),
			m_scc_rank_in(),
			m_weak_coloring(),
			m_weakly_connected_components(),
			m_heuristic_out(),
			m_heuristic_in(),
			m_basic_topological_ordering(),
			m_basic_topological_ordering_in(),
			m_basic_topological_ordering_inverse(),
			m_basic_topological_ordering_inverse_in(),
			m_scc_big_components(),
			m_node_names(vnames),
			m_namemap()
{
	for (node_t i = 0; i < vnames.size(); ++i)
    {
        m_namemap[vnames[i]] = i;
    }
    
}

void DiGraph::add_edge(node_t from, node_t to, weight_t weight)
{
// 	if (m_edge_values(from,to) == 0)
// 	{
		m_outgraph[from].emplace_back(to,weight);
		m_ingraph[to].emplace_back(from,weight);
// 	}
// 	m_edge_values(from,to) = weight;
	m_processed = false;
}

void DiGraph::add_edge(const string& from, const string& to, weight_t weight)
{
    add_edge(m_namemap[from],m_namemap[to],weight);
}


void DiGraph::process()
{
	if (!m_processed)
	{
		cout << "Processing graph..." << endl;
		Chronometer C;
		find_weakly_connected_components();
		cout << "weak coloring = " << m_weak_coloring << endl;
		cout << "before removing, size: " << num_vertices() << endl;
		remove_bad_nodes();
		cout << "After removing, new size: " << m_n << endl;
		
		find_strongly_connected_components();
		
		cout << "Scc: " << m_scc_coloring << endl;
		
		heuristic_processing();
		
		cout << "Heuristic: " << m_heuristic_out << endl;
		m_processed = true;
		cout << "Finished processing graph after " << C.Peek() << " seconds" << endl;
	}
}

Path DiGraph::forward_backward_dfs(node_t node, double mnswi)
{
	Path FP = dfs_search_path_forward(node,mnswi);
	
	Path BP = FP;
	
	int numremoves  = std::min(20,int(BP.size()/2)-1);
	for (int i = 0; i < numremoves; ++i) BP.pop_front();
	
	dfs_search_path_reverse(BP,mnswi);
    if (BP.value() > FP.value())
        FP = BP;
    
    return FP;
}

Path DiGraph::dfs_search(double mnswi, int numrestarts) const
{
	Path Best;
	
	if (numrestarts > num_vertices())
		numrestarts = num_vertices();
	
    for (int i = 0; i < numrestarts; ++i)
    {
		node_t node = m_basic_topological_ordering[i/2];
		if (i%2 == 1)
			node = rand()%num_vertices();
		
        Path P = dfs_search_path_forward(node,mnswi);
        if (P.value() > Best.value())
            Best = std::move(P);
    }
    
    return Best;
}


Path DiGraph::FindLongestSimplePath(double numseconds, double numsecondsDFSnoimprovement, int numrestarts)
{
    Chronometer C;
	process();
	cout << "Doing DFS search..." << endl;
	Path best = dfs_search(numsecondsDFSnoimprovement,numrestarts);
	double timeleft = numseconds-C.Peek();
	
	cout << "Done DFS search! Value = " << best.value() << endl;
	cout << "Time taken for dfs: " << C.Peek() << endl;
	
	cout << "Doing PTO improving search for " << timeleft << "s" << endl;
	pto_search(best,timeleft);
	
	return best;
}

void DiGraph::pto_search(Path& A, double maxnumseconds) const
{
	Chronometer C;
	
    PseudoTopoOrder PTO = get_random_pseudotopological_order();
// 	cout << "Before applying A, PTO.Value() = " << PTO.Value() << endl;
	
	PTO.apply(A);
	
// 	cout << "After applying A, PTO.Value() = " << PTO.Value() << endl;

	PTO.open_edges_until_no_more_improvement_found(maxnumseconds);
	
	cout << "After opening edges, PTO.Value() = " << PTO.Value() << endl;
	
    A = PTO.get_path();
}

size_t DiGraph::num_edges() const
{
	size_t toReturn = 0;
	for (node_t i = 0; i < m_n; ++i)
	{
		toReturn += m_outgraph[i].size();
	}
	return toReturn;
}


void DiGraph::remove_bad_nodes()
{
	vector<node_t> toRemove;
	toRemove.reserve(m_n/2);
	
// 		cout << "removing bad nodes!" << endl;
// 	int val = dfs_search_path(0.1).Value();
	
	size_t val = 0;
	for (size_t i = 0; i < m_weakly_connected_components.size(); ++i)
	{
		if (val < m_weakly_connected_components[i].size())
			val = m_weakly_connected_components[i].size();
	}

	for (size_t i = 0; i < m_weakly_connected_components.size(); ++i)
	{
		if (m_weakly_connected_components[i].size() < val)
		{	
			for (auto x : m_weakly_connected_components[i])
				toRemove.push_back(x);
		}
	}
// 	cout << "done adding toremove! sorting..." << endl;
// 	sort(toRemove.begin(), toRemove.end());
	
// 	for (int i = 0; i < n; ++i)
// 	{
// 		if (outgraph[i].size() == 0 && ingraph[i].size() == 0)
// 			toRemove.push_back(i);
// 	}
// 	cout << "Removing: " << toRemove.size() << " nodes" << endl;
	remove_nodes(toRemove);
}


DiGraph DiGraph::with_nodes_removed(vector<node_t>& toRemove) const
{
    if (toRemove.empty())
	{
		return *this;
	}
	sort(toRemove.begin(), toRemove.end());
	
	size_t new_n = m_n-toRemove.size();
// 	DiGraph D(new_n);
	vector<string> new_names;
	new_names.reserve(new_n);
    
	vector<node_t> removalfunction(m_n,INVALID_NODE);
	vector<node_t> removalfunctioninverse(new_n,0);
	
	node_t j = 0;
	node_t w = 0;
	for (node_t i = 0; i < m_n; ++i)
	{
		if (i != toRemove[w])
		{
			removalfunction[i] = j;
			removalfunctioninverse[j] = i;
			new_names.emplace_back(m_node_names[i]);
			++j;
		} else
		{
			++w;
		}
	}
	
	DiGraph D(new_names);
	
	for (node_t v = 0; v < new_n; ++v)
    {
        auto oldv = removalfunctioninverse[v];
        for (auto oldneigh : outneighbors(oldv))
        {
            auto newneigh = removalfunction[oldneigh];
            D.add_edge(v,newneigh,oldneigh.Weight());
        }
    }
    return D;
}

void DiGraph::remove_nodes(vector<node_t>& toRemove)
{
	*this = with_nodes_removed(toRemove);
}

void DiGraph::dfs_search_path_forward(Path& P, double maxnumseconds) const
{
	ExpandGreedyBack(*this,P);
// 	cout << "(" << ChronometerPeek() << ", " << P.Value() << ")," << endl;
	
	Chronometer C;
	Path Q = P;
	auto comp = [this](node_t a, node_t b)
	{
		return ex_compare(a,b);
	};
	while (C.Peek() < maxnumseconds && dfs_outnext(*this,Q,comp))
	{
		if (Q.value() > P.value())
		{
			P = Q;
			C.Reset();
// 			cout << "(" << C.Peek() << ", " << P.value() << ")," << endl;
// 			cout << P.value()/2 << endl;
		}
	}
}

void DiGraph::dfs_search_path_reverse(Path& P, double maxnumseconds) const
{
	ExpandGreedyBack(*this,P);
// 	cout << "(" << ChronometerPeek() << ", " << P.Value() << ")," << endl;
	
	Chronometer C;
	Path Q = P;
	auto comp = [this](node_t a, node_t b)
	{
		return in_compare(a,b);
	};
	while (C.Peek() < maxnumseconds && dfs_innext(*this,Q,comp))
	{
		if (Q.value() > P.value())
		{
			P = Q;
			C.Reset();
// 			cout << "(" << ChronometerPeek() << ", " << P.Value() << ")," << endl;
// 			cout << P.value()/2 << endl;
		}
	}
}

Path DiGraph::dfs_search_path_forward(node_t start, double maxnumseconds) const
{
	Path P(num_vertices(),start);
	dfs_search_path_forward(P,maxnumseconds);
	return P;
}

Path DiGraph::dfs_search_path_reverse(node_t start, double maxnumseconds) const
{
	Path P(num_vertices(),start);
	dfs_search_path_reverse(P,maxnumseconds);
	return P;
}

void DiGraph::heuristic_processing()
{
	m_heuristic_out.resize(m_n,0);
	m_heuristic_in.resize(m_n,0);
	
	// 	double maxtime = 0.2/n;
	
// 	#pragma omp parallel for
	for (node_t i = 0; i < m_n; ++i)
	{
		m_heuristic_out[i] = get_heuristic_out(i);
		m_heuristic_in[i] = get_heuristic_in(i);
	}

	
	m_basic_topological_ordering = range<node_t>(m_n);
	sort (m_basic_topological_ordering.begin(), m_basic_topological_ordering.end(), [this](node_t a, node_t b) -> bool
	{
		//First, order by rank
		if (rank_out(a) < rank_out(b))
			return false;
		if (rank_out(a) > rank_out(b))
			return true;
		
		//If someone doesn't have out-neighbours, they should go last
		if (outdegree(a) == 0)
			return false;
		if (outdegree(b) == 0)
			return true;
		
		//The ones who have NO OTHER entry points should go first
		if (indegree(b) == 1)
			return false;
		if (indegree(a) == 1 )
			return true;
		
		//If all else fails, sort by heuristic
		return m_heuristic_out[a] < m_heuristic_out[b];
	});
	
	m_basic_topological_ordering_inverse.resize(m_n);
	for (node_t i = 0; i < m_n; ++i)
	{
		m_basic_topological_ordering_inverse[m_basic_topological_ordering[i]] = i; 
	}
	
	m_basic_topological_ordering_in = range<node_t>(m_n);
	sort (m_basic_topological_ordering_in.begin(), m_basic_topological_ordering_in.end(), [this](node_t a, node_t b) -> bool
	{
		if (rank_in(a) > rank_in(b))
			return false;
		if (rank_in(a) < rank_in(b))
			return true;
		
		
		if (indegree(a) == 0)
			return false;
		if (indegree(b) == 0)
			return true;
		
		if (outdegree(b) == 1)
			return false;
		if (outdegree(a) == 1 )
			return true;

		return m_heuristic_in[a] < m_heuristic_in[b];
	});

	m_basic_topological_ordering_inverse_in.resize(m_n);
	for (node_t i = 0; i < m_n; ++i)
	{
		m_basic_topological_ordering_inverse_in[m_basic_topological_ordering_in[i]] = i; 
	}
	
	
	for (size_t i = 0; i < m_n; ++i)
	{
		sort(m_outgraph[i].begin(), m_outgraph[i].end(), [this] (node_t a, node_t b) -> bool
		{
			return ex_compare(a,b);
		});
		
		sort(m_ingraph[i].begin(), m_ingraph[i].end(), [this] (node_t a, node_t b) -> bool
		{
			return in_compare(a,b);
		});
	}
	
}

void DiGraph::DFSUtil(node_t v, vector<bool>& visited)
{
    visited[v] = true;
 
    for (auto i : m_outgraph[v])
	{
        if (!visited[i])
            DFSUtil(i, visited);
	}
}
 
 
 void DiGraph::DFSUtilReversed(node_t v, vector< char >& visited, int current)
{
    visited[v] = true;
// 	cout << " v = " << v << " and current = " << current << " and scc.size() = " << strongly_connected_components.size() << endl;
	m_strongly_connected_components[current].push_back(v);
	m_scc_coloring[v] = current;
// 	cout << "good" << endl;
    for (auto i : m_ingraph[v])
	{
        if (!visited[i])
            DFSUtilReversed(i, visited,current);
	}
}
 
 
void DiGraph::topo_fill_order(node_t v, vector< char >& visited, stack< node_t >& Stack)
{
    // Mark the current node as visited and print it
    visited[v] = true;
 
    // Recur for all the vertices outgraphacent to this vertex
    for(auto i : m_outgraph[v])
        if(!visited[i])
            topo_fill_order(i, visited, Stack);
 
    // All vertices reachable from v are processed by now, push v 
    Stack.push(v);
}

void DiGraph::find_strongly_connected_components()
{
	m_strongly_connected_components.clear();
	m_strongly_connected_components.reserve(m_n/4);

	m_scc_coloring.resize(m_n);
	
    stack<node_t> Stack;
 
    vector<char> visited(m_n,0);
 
    for(int i = 0; i < m_n; ++i)
        if(visited[i] == false)
            topo_fill_order(i, visited, Stack);
 
    for(int i = 0; i < m_n; i++)
        visited[i] = false;
	
	int current = 0;
// 	cout << "before starting, scc.size() = " << m_strongly_connected_components.size() << endl;
    while (Stack.empty() == false)
    {
        int v = Stack.top();
        Stack.pop();
 
		m_strongly_connected_components.push_back({});
        if (visited[v] == false)
        {
            DFSUtilReversed(v, visited,current);
// 			cout << endl;
            ++current;
			m_strongly_connected_components.push_back({});
        }
    }
    
//     cout << "Finished first loop" << endl;
	
    while (m_strongly_connected_components.back().empty())
		m_strongly_connected_components.pop_back();
	
//     cout << "Finished popping the empty components. WTF?" << endl;

	int num_scc = m_strongly_connected_components.size();
	m_scc_rank_in.resize(m_n,0);
	
	for (int scc = 0; scc < num_scc; ++scc)
	{
		for (int x : m_strongly_connected_components[scc])
		{
			// x is the actual node in the scc component
			for (auto v : m_ingraph[x])
			{
				//v is the actual node in the scc component
				int connected_component_neigh = m_scc_coloring[v];
				if (connected_component_neigh == scc)
					continue;
				int candidate = m_scc_rank_in[connected_component_neigh] + 1;
				if (candidate > m_scc_rank_in[scc])
					m_scc_rank_in[scc] = candidate;
			}
		}
	}
	
// 	cout << "Finished with this weird loop." << endl;
	
	m_scc_rank_out.resize(num_scc,0);

	for (int scc = num_scc-1; scc >= 0; --scc)
	{
		for (int x : m_strongly_connected_components[scc])
		{
			// x is the actual node in the scc component
			for (auto v : m_outgraph[x])
			{
				//v is the actual node in the scc component
				int connected_component_neigh = m_scc_coloring[v];
				if (connected_component_neigh == scc)
					continue;
				int candidate = m_scc_rank_out[connected_component_neigh] + 1;
				if (candidate > m_scc_rank_out[scc])
					m_scc_rank_out[scc] = candidate;
			}
		}
	}
// 	cout << "Finished with this other weird loop." << endl;

	int start = 0;
	for (int color = 0; color < num_scc; ++color)
	{
		int s = m_strongly_connected_components[color].size();
		if (s > 1)
		{
			m_scc_big_components.push_back(SccInfo(start,start+s,color));
		}
		start += s;
	}

// 	cout << "Finished with this third weird loop." << endl;

// 	for (int i = 0; i < m_scc_big_components.size(); ++i)
// 	{
// 		auto u = m_scc_big_components[i];
// 		cout << i << ": " << u.start << " --> " << u.end << " w color " << u.color << endl;
// 	}
	
// 	cout << "Finished printing big component info." << endl;
// 	cout << " with m_strongly_connected_components.size() = " << m_strongly_connected_components.size() << endl;
// 	cout << " and num_scc = " << num_scc << endl;
// 	for (int i = 0; i < num_scc; ++i)
// 	{
// 		cout << i << endl;
// 		auto x = m_strongly_connected_components[i];
// 		if (x.size() > 1)
// 			cout << i << ": " << x << endl;
// 	}
	
// 	cout << "out of..." << m_strongly_connected_components.size();
// 	cout << "scc_rank_out = " << m_scc_rank_out << endl;
	
}

PseudoTopoOrder DiGraph::get_random_pseudotopological_order() const
{
	vector<node_t> topo_sort(m_n,0); // at index i, it has a node topo_sort[i]
	vector<node_t> topo_sort_inverse(m_n,0); // at node i, x=topo_sort_inverse[i] is such that topo_sort[x] = i
	int i = 0;
	for (auto X : m_strongly_connected_components)
	{
		random_shuffle(X.begin(), X.end());
		for (auto x : X)
		{
			topo_sort[i] = x;
			topo_sort_inverse[x] = i;
			++i;
		}
	}
	return PseudoTopoOrder(*this, std::move(topo_sort),std::move(topo_sort_inverse));
}

void DiGraph::DFSUtilWeak(node_t node, int color)
{
	m_weak_coloring[node] = color;
	m_weakly_connected_components.back().push_back(node);
// 	m_weakly_connected_components_values.back() += vertex_values[node]-1;
	for (auto x : m_outgraph[node])
	{
		if (m_weak_coloring[x] != -1)
			continue;
		DFSUtilWeak(x,color);
	}
	
	for (auto x : m_ingraph[node])
	{
		if (m_weak_coloring[x] != -1)
			continue;
		DFSUtilWeak(x,color);
	}
}

void DiGraph::find_weakly_connected_components()
{
	m_weak_coloring.clear();
	m_weak_coloring.resize(m_n,-1);
	m_weakly_connected_components.clear();
// 	m_weakly_connected_components_values.clear();
	int minvalidcoloring = 0;
	for (int start = 0; start < m_n; ++start)
	{
		if (m_weak_coloring[start] != -1)
			continue;
		m_weakly_connected_components.push_back({});
// 		m_weakly_connected_components_values.push_back(0);
		DFSUtilWeak(start, minvalidcoloring);
		
		++minvalidcoloring;
	}
	
}


double DiGraph::get_heuristic_out(node_t node)
{
    double a1 = m_params[0];
	double a2 = m_params[1];
	double a3 = m_params[2];
	double a4 = m_params[3];
	double heuristicex = 0;
	
	for (auto x : m_outgraph[node])
	{
		heuristicex += a1;
		for (auto y : m_outgraph[x])
		{
			heuristicex += a2;
			for (auto z : m_outgraph[y])
			{
				heuristicex += a3+a4*m_outgraph[z].size();
// 				for (auto r : outgraph[z])
// 				{
// 					heuristicex += a4+a5*outgraph[r]./*size()*/;
// 				}
			}
		}
	}
	return heuristicex;
}


double DiGraph::get_heuristic_in(node_t node)
{
// 	return 0;
	double a1 = m_params[4];
	double a2 = m_params[5];
	double a3 = m_params[6];
	double a4 = m_params[7];
	double heuristicin = 0;

	for (auto x : m_ingraph[node])
	{
		heuristicin += a1;
		for (auto y : m_ingraph[x])
		{
			heuristicin  += a2+a3*m_ingraph[y].size();
			for (auto z : m_ingraph[y])
			{
				heuristicin  += a3+a4*m_ingraph[z].size(); 
			}
		}
	}
	return heuristicin;
}

DiGraph DiGraph::CreateRandomDiGraph(int n, double p)
{
	DiGraph D(n);
	
	for (int i = 0; i < n; ++i)
	{
		for (int j = i+1; j < n; ++j)
		{
			if (probability_of_true(p))
			{
				int a = i;
				int b = j;
				if (rand()%2 == 1) swap(a,b);
				D.add_edge(a,b);
			}
		}
	}
// 	D.process();
	return D;
}


DiGraph DiGraph::CreateRandomWeightedDiGraph(int n, double p, weight_t minweight, weight_t maxweight)
{
	DiGraph D(n);
	
	for (int i = 0; i < n; ++i)
	{
		for (int j = i+1; j < n; ++j)
		{
			if (probability_of_true(p))
			{
				int a = i;
				int b = j;
				weight_t w = random_real(minweight,maxweight);
				if (rand()%2 == 1) swap(a,b);
				D.add_edge(a,b,w);
			}
		}
	}
// 	D.process();
	return D;
}

std::ostream& operator<<(std::ostream& os, const ParamType& a)
{
	os << endl << "\tex: ";
	for (int i = 0; i < 4; ++i)
		os << a[i] << " ";
	os << endl;
	
	os << "\tin: ";
	for (int i = 4; i < 8; ++i)
		os << a[i] << " ";
	os << endl;
	return os;
}

std::ostream& operator<<(std::ostream& os, const DiGraph& M)
{
    os << "Digraph on " << M.num_vertices() << " vertices: " << M.get_vertex_names();
    for (int i = 0; i < M.num_vertices(); ++i)
    {
        string name = M.get_vertex_name(i);
        for (auto v : M.outneighbors(i))
        {
			os << endl;
            string vname = M.get_vertex_name(v);
            os << name << u8" ⟼ " << vname << " with weight " << double(v.Weight());
        }
    }
    return os;
}

void ExpandGreedyBack(const DiGraph& G, Path& P)
{
	while (true)
	{
		auto l = P.back();
		auto& Neighs = G.outneighbors(l);
		
		auto t = P.first_not_explored(Neighs);
		if (t == INVALID_NODE)
			break;
		P.emplace_back(t,t.Weight());
	}
}


void ExpandGreedyFront(const DiGraph& G, Path& P)
{
	while (true)
	{
		auto l = P.front();
		auto& Neighs = G.inneighbors(l);
		
		auto t = P.first_not_explored(Neighs);
		if (t == INVALID_NODE)
			break;
		P.emplace_front(t,t.Weight());
	}
}
