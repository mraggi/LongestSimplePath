#include "digraph.hpp"
#include "path.hpp"
#include "pseudotopoorder.hpp"
#include <sstream>

//Heuristics

DiGraph::DiGraph(node_t numNodes) :
			m_n(numNodes),
			m_exgraph(numNodes),
			m_ingraph(numNodes),
			m_edge_values(numNodes),
			m_processed(false),
			m_strongly_connected_components(),
			m_scc_coloring(),
			m_scc_rank_ex(),
			m_scc_rank_in(),
			m_weak_coloring(),
			m_weakly_connected_components(),
			m_heuristic_ex(),
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
			m_exgraph(vnames.size()),
			m_ingraph(vnames.size()),
			m_edge_values(vnames.size()),
			m_processed(false),
			m_strongly_connected_components(),
			m_scc_coloring(),
			m_scc_rank_ex(),
			m_scc_rank_in(),
			m_weak_coloring(),
			m_weakly_connected_components(),
			m_heuristic_ex(),
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
	if (m_edge_values(from,to) == 0)
	{
		m_exgraph[from].push_back(to);
		m_ingraph[to].push_back(from);
	}
	m_edge_values(from,to) = weight;
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
		clock_t start = clock();
// 		cout << "Processing graph..." << endl;
		find_weakly_connected_components();
// 		cout << "weak coloring = " << m_weak_coloring << endl;
		remove_bad_nodes();
// 		cout << "After removing, new size: " << m_n << endl;
		
		find_strongly_connected_components();
		
// 		cout << "Scc: " << m_scc_coloring << endl;
		
		heuristic_processing();
		
// 		cout << "Heuristic: " << m_heuristic_ex << endl;
		m_processed = true;
// 		cout << "Finished processing after " << diffclock(clock(),start) << " seconds" << endl;
	}
}

Path DiGraph::dfs_search(double maxnumsecondswithoutimprovement) const
{
    clock_t start = clock();
    Path A = dfs_search_path_forward(m_basic_topological_ordering[0],maxnumsecondswithoutimprovement);
    node_t i = 1;
    double totaltime = 5*maxnumsecondswithoutimprovement;
    double timeleft = totaltime - diffclock(clock(),start);
    while (i < num_vertices() && timeleft > maxnumsecondswithoutimprovement)
    {
        Path P = dfs_search_path_forward(m_basic_topological_ordering[i],maxnumsecondswithoutimprovement);
        if (P.Value() > A.Value())
            A = P;
        ++i;
        timeleft = totaltime - diffclock(clock(),start);
    }
    Path P = A;
    for (int i = 0; i < 15; ++i) A.PopFront();
	dfs_search_path_reverse(A,maxnumsecondswithoutimprovement/3.0);
    if (P.Value() > A.Value())
        A = P;
    
    return A;
}


Path DiGraph::FindLongestSimplePath(double numseconds)
{
	process(); // do this BEFORE starting the count.
    
    Chronometer();
	clock_t start = clock();
	
	double dfstime = 0.07 + numseconds/50.0;
	
	Path best = dfs_search(dfstime);
	
	double elapsed = diffclock(clock(),start);
	double timeleft = numseconds-elapsed;
	
// 	cout << "Done DFS search! Value = " << best.Value() << endl;
	
	pto_search(best,timeleft);
	
	return best;
}

void DiGraph::pto_search(Path& A, double maxnumseconds) const
{
   	clock_t start = clock();

    PseudoTopoOrder PTO = get_random_pseudotopological_order();
// 	cout << "Before applying A, PTO.Value() = " << PTO.Value() << endl;
	
	PTO.apply(A);
	
// 	cout << "After applying A, PTO.Value() = " << PTO.Value() << endl;
// 	cout << "(" << ChronometerPeek() << ", " << PTO.Value() << ")," << endl;

	double timeleft = maxnumseconds-diffclock(clock(),start);
// 	cout << "*****************After second DFS, I have " << timeleft << "s of time left" << endl;
	
	PTO.open_edges_until_no_more_improvement_found(min(timeleft,5.0));
	
// 	cout << "After opening edges, PTO.Value() = " << PTO.Value() << endl;
	
    if (A.Value() < PTO.Value())
        A = PTO.get_path();
	
// 	cout << "which should be the same as the value of A: " << A.Value() << endl;
	
	timeleft = maxnumseconds-diffclock(clock(),start);

// 	cout << "After opening edges, I have " << timeleft << "s of time left" << endl;
	Path WA = A;
	while (timeleft > 0.1)
	{
		int numtopop = rand()%8+8;
		for (int i = 0; i < numtopop; ++i)
		{
			WA.PopBack();
			WA.PopFront();
		}
		PTO.randomize();
		PTO.apply(WA);
		PTO.open_edges_until_no_more_improvement_found(min(timeleft,1.5));
		if (PTO.Value() > A.Value())
		{
// 			cout << "holy shit, I advanced here to " << PTO.Value() << endl;
			A = PTO.get_path();
		}
		WA = A;
        timeleft = maxnumseconds-diffclock(clock(),start);
	}
}

Path DiGraph::FindLongestSimplePathPureDFS(double numseconds)
{
	process();
	
	Path best(this,0);
	
	clock_t start = clock();
	Chronometer();
	
	double dfstime = numseconds/2.05 - 0.1;
	
	Path A = dfs_search_path_forward(m_basic_topological_ordering[0],dfstime);
// 	cout << "Done forward search! Value = " << A.Value() << endl;
	
	double elapsed = diffclock(clock(),start);
	double timeleft = numseconds-elapsed;
// 	cout << "After the first DFS, I have " << timeleft << "s of time left" << endl;
	
	if (timeleft > dfstime)
	{
		for (int i = 0; i < max(min(12,int(timeleft*20)),1); ++i) A.PopFront();
		
// 		cout << "Starting backward search" << endl;
		dfs_search_path_reverse(A,dfstime);
	}
	
	if (A.Value() > best.Value())
		best = A;
	
	elapsed = diffclock(clock(),start);
	timeleft = numseconds-elapsed;
// 	cout << "After opening edges, I have " << timeleft << "s of time left" << endl;
	
	return best;
}

size_t DiGraph::num_edges() const
{
	size_t toReturn = 0;
	for (node_t i = 0; i < m_n; ++i)
	{
		toReturn += m_exgraph[i].size();
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
	sort(toRemove.begin(), toRemove.end());
	
// 	for (int i = 0; i < n; ++i)
// 	{
// 		if (exgraph[i].size() == 0 && ingraph[i].size() == 0)
// 			toRemove.push_back(i);
// 	}
// 	cout << "Removing: " << toRemove.size() << " nodes" << endl;
	remove_nodes(toRemove);
}

void DiGraph::remove_nodes(const vector<node_t>& toRemove)
{
	if (toRemove.empty())
	{
		return;
	}
	
	size_t new_n = m_n-toRemove.size();
	
	vector<node_t> removalfunction(m_n,INVALID_NODE);
	vector<node_t> removalfunctioninverse(new_n,0);
	
	vector<string> new_movienames(new_n);
	vector<vector<node_t>> new_exgraph(new_n);
	vector<vector<node_t>> new_ingraph(new_n);
	SquareSparseMatrix<weight_t> new_edge_values(new_n);
    
    vector<string> new_names(new_n);
    unordered_map<string,node_t> new_namemap;
    
	int j = 0;
	int w = 0;
	for (node_t i = 0; i < m_n; ++i)
	{
		if (i != toRemove[w])
		{
			removalfunction[i] = j;
			removalfunctioninverse[j] = i;
			++j;
		} else
		{
			++w;
		}
	}
	
	for (node_t j = 0; j < new_n; ++j)
	{
		auto i = removalfunctioninverse[j];
        new_names[j] = m_node_names[i];
        new_namemap[new_names[j]] = j;
		for (auto x : m_exgraph[i])
		{
			auto new_x = removalfunction[x];
			new_exgraph[j].push_back(new_x);
			new_ingraph[new_x].push_back(j);
			new_edge_values(j,new_x) = m_edge_values(i,x);
		}
	}
// 	cout << "Finished setting exgraph and ingraph. Copying..." << endl;
	m_n = new_n;
	m_exgraph = std::move(new_exgraph);
	m_ingraph = std::move(new_ingraph);
	m_edge_values = std::move(new_edge_values);
    m_node_names = std::move(new_names);
    m_namemap = std::move(new_namemap);
}

void DiGraph::dfs_search_path_forward(Path& P, double maxnumseconds) const
{
	P.ExpandGreedyBack();
// 	cout << "(" << ChronometerPeek() << ", " << P.Value() << ")," << endl;
	
	clock_t start = clock();
	Path Q = P;
	while (diffclock(clock(),start) < maxnumseconds && Q.DFSBackNext())
	{
		if (Q.Value() > P.Value())
		{
			P = Q;
			start = clock();
// 			cout << "(" << ChronometerPeek() << ", " << P.Value() << ")," << endl;
		}
	}
}

void DiGraph::dfs_search_path_reverse(Path& P, double maxnumseconds) const
{
	P.ExpandGreedyFront();
// 	cout << "(" << ChronometerPeek() << ", " << P.Value() << ")," << endl;
	
	clock_t start = clock();
	Path Q = P;
	while (diffclock(clock(),start) < maxnumseconds && Q.DFSFrontNext())
	{
		if (Q.Value() > P.Value())
		{
			P = Q;
			start = clock();
// 			cout << "(" << ChronometerPeek() << ", " << P.Value() << ")," << endl;

		}
	}
}

Path DiGraph::dfs_search_path_forward(node_t start, double maxnumseconds) const
{
	Path P(this,start);
	dfs_search_path_forward(P,maxnumseconds);
	return P;
}

Path DiGraph::dfs_search_path_reverse(node_t start, double maxnumseconds) const
{
	Path P(this,start);
	dfs_search_path_reverse(P,maxnumseconds);
	return P;
}

void DiGraph::heuristic_processing()
{
	m_heuristic_ex.resize(m_n,0);
	m_heuristic_in.resize(m_n,0);
	
	// 	double maxtime = 0.2/n;
	
// 	#pragma omp parallel for
	for (node_t i = 0; i < m_n; ++i)
	{
		m_heuristic_ex[i] = get_heuristic_ex(i);
		m_heuristic_in[i] = get_heuristic_in(i);
	}

	
	m_basic_topological_ordering = range<node_t>(m_n);
	sort (m_basic_topological_ordering.begin(), m_basic_topological_ordering.end(), [this](node_t a, node_t b) -> bool
	{
		if (rank_ex(a) > rank_ex(b))
			return true;
		if (rank_ex(a) < rank_ex(b))
			return false;
		
// 		if (m_heuristic_ex[a] == 0 && m_heuristic_ex[b] == 0)
// 		{
// 			return m_vertex_values[a] > m_vertex_values[b];
// 		}
		
		if (m_heuristic_ex[a] == 0)
			return false;
		if (m_heuristic_ex[b] == 0)
			return true;
		
// 		if (ingraph[a].size() == 1 && ingraph[b].size() == 1)
// 			return heuristic_ex[a] < heuristic_ex[b];
		
		if (m_ingraph[a].size() == 1 )
			return true;
		if (m_ingraph[b].size() == 1)
			return false;
		
		return m_heuristic_ex[a] < m_heuristic_ex[b];
	});
// 	random_shuffle(basic_topological_ordering.begin(), basic_topological_ordering.end());
	m_basic_topological_ordering_inverse.resize(m_n);
	for (node_t i = 0; i < m_n; ++i)
	{
		m_basic_topological_ordering_inverse[m_basic_topological_ordering[i]] = i; 
	}
	
	m_basic_topological_ordering_in = range<node_t>(m_n);
	sort (m_basic_topological_ordering_in.begin(), m_basic_topological_ordering_in.end(), [this](node_t a, node_t b) -> bool
	{
		if (rank_in(a) < rank_in(b))
			return true;
		if (rank_in(a) > rank_in(b))
			return false;
		
// 		if (m_heuristic_in[a] == 0 && m_heuristic_in[b] == 0)
// 		{
// 			return vertex_values[a] > vertex_values[b];
// 		}
		
		if (m_heuristic_in[a] == 0)
			return false;
		if (m_heuristic_in[b] == 0)
			return true;
		
// 		if (ingraph[a].size() == 1 && ingraph[b].size() == 1)
// 			return heuristic_ex[a] < heuristic_ex[b];
		
		if (m_exgraph[a].size() == 1 )
			return true;
		if (m_exgraph[b].size() == 1)
			return false;
		
		return m_heuristic_in[a] < m_heuristic_in[b];
	});
// 	random_shuffle(basic_topological_ordering.begin(), basic_topological_ordering.end());
	m_basic_topological_ordering_inverse_in.resize(m_n);
	for (node_t i = 0; i < m_n; ++i)
	{
		m_basic_topological_ordering_inverse_in[m_basic_topological_ordering_in[i]] = i; 
	}
	
	
	for (size_t i = 0; i < m_n; ++i)
	{
// 		random_shuffle(exgraph[i].begin(), exgraph[i].end());
// 		random_shuffle(ingraph[i].begin(), ingraph[i].end());
		sort(m_exgraph[i].begin(), m_exgraph[i].end(), [this] (node_t a, node_t b) -> bool
		{
			return ex_compare(a,b);
		});
		
		sort(m_ingraph[i].begin(), m_ingraph[i].end(), [this] (node_t a, node_t b) -> bool
		{
			return in_compare(a,b);
		});
	}
	
}

void DiGraph::DFSUtil(int v, vector<bool>& visited)
{
    visited[v] = true;
 
    for (auto i : m_exgraph[v])
	{
        if (!visited[i])
            DFSUtil(i, visited);
	}
}
 
 
 void DiGraph::DFSUtilReversed(int v, vector< char >& visited, int current)
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
 
 
void DiGraph::topo_fill_order(int v, vector< char >& visited, stack< node_t >& Stack)
{
    // Mark the current node as visited and print it
    visited[v] = true;
 
    // Recur for all the vertices exgraphacent to this vertex
    for(auto i : m_exgraph[v])
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
	
	m_scc_rank_ex.resize(num_scc,0);

	for (int scc = num_scc-1; scc >= 0; --scc)
	{
		for (int x : m_strongly_connected_components[scc])
		{
			// x is the actual node in the scc component
			for (auto v : m_exgraph[x])
			{
				//v is the actual node in the scc component
				int connected_component_neigh = m_scc_coloring[v];
				if (connected_component_neigh == scc)
					continue;
				int candidate = m_scc_rank_ex[connected_component_neigh] + 1;
				if (candidate > m_scc_rank_ex[scc])
					m_scc_rank_ex[scc] = candidate;
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
// 	cout << "scc_rank_ex = " << m_scc_rank_ex << endl;
	
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

void DiGraph::DFSUtilWeak(int node, int color)
{
	m_weak_coloring[node] = color;
	m_weakly_connected_components.back().push_back(node);
// 	m_weakly_connected_components_values.back() += vertex_values[node]-1;
	for (auto x : m_exgraph[node])
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


double DiGraph::get_heuristic_ex(int node)
{
    double a1 = m_params[0];
	double a2 = m_params[1];
	double a3 = m_params[2];
	double a4 = m_params[3];
	double heuristicex = 0;
	
	for (auto x : m_exgraph[node])
	{
		heuristicex += a1;
		for (auto y : m_exgraph[x])
		{
			heuristicex += a2;
			for (auto z : m_exgraph[y])
			{
				heuristicex += a3+a4*m_exgraph[z].size();
// 				for (auto r : exgraph[z])
// 				{
// 					heuristicex += a4+a5*exgraph[r]./*size()*/;
// 				}
			}
		}
	}
	return heuristicex;
}


double DiGraph::get_heuristic_in(int node)
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

Path DiGraph::get_random_path(double maxnumseconds) const
{
	node_t startnode = rand()%m_n;
	Path trivial(this,startnode);
	Path toReturn = trivial;

	vector<Path> frontier = {trivial};
	clock_t st = clock();
	while (!frontier.empty() && diffclock(clock(), st) < 0.5*maxnumseconds)
	{
		Path P = frontier.back();
		frontier.pop_back();
		if (P.Value() > toReturn.Value())
			toReturn = P;
		
		int lastnode = P.get_path().back();
		auto exx = m_exgraph[lastnode];
		random_shuffle(exx.begin(), exx.end());
		sort(exx.begin(), exx.end(), [this](int a, int b) -> bool{
			return m_scc_rank_ex[a] < m_scc_rank_ex[b];
		});
		for (auto x : exx)
		{
			if (P.IsNodeInPath(x) == false)
			{
				
				Path PP = P;
				PP.AddNodeBack(x);
				frontier.push_back(PP);
				
			}
		}
	}
// 	if (frontier.empty())
// 	{
// 		return get_random_path()(maxnumseconds-diffclock(clock(), st));
// 	}
	vector<Path> frontier2;
	frontier2.push_back(toReturn);
	while (!frontier2.empty() && diffclock(clock(), st) < maxnumseconds)
	{
		Path P = frontier2.back();
		frontier2.pop_back();
		if (P.Value() > toReturn.Value())
			toReturn = P;
		int lastnode = P.get_path().front();
		auto inn = m_ingraph[lastnode];
		random_shuffle(inn.begin(), inn.end());
		sort(inn.begin(), inn.end(), [this](int a, int b) -> bool{
			return m_scc_rank_in[a] < m_scc_rank_in[b];
		});
		for (auto x : inn)
		{
			if (P.IsNodeInPath(x) == false)
			{
				
				Path PP = P;
				PP.AddNodeFront(x);
				frontier2.push_back(PP);
			}
		}
	}
	
	return toReturn;
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
    cout << "Digraph on " << M.num_vertices() << " vertices: " << M.get_vertex_names() << endl;
    for (int i = 0; i < M.num_vertices(); ++i)
    {
        string name = M.get_vertex_name(i);
        for (auto v : M.exneighbors(i))
        {
            string vname = M.get_vertex_name(v);
            os << name << u8" ⟼ " << vname << " with weight " << double(M.edge_value(i,v)) << endl;
        }
    }
    return os;
}

