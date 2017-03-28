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
		m_node_names[i] = std::to_string(i);
		m_namemap[std::to_string(i)] = i;
	}
	

}

DiGraph::DiGraph(const std::vector<std::string>& vertex_names)  :
	m_n(vertex_names.size()),
	m_outgraph(m_n),
	m_ingraph(m_n),
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
	m_node_names(vertex_names),
	m_namemap()
{
	for (node_t i = 0; i < m_node_names.size(); ++i)
	{
		m_namemap[m_node_names[i]] = i;
	}

}

void DiGraph::add_edge(node_t from, node_t to, weight_t weight)
{
// 	if (m_edge_values(from,to) == 0)
// 	{
	m_outgraph[from].emplace_back(to, weight);
	m_ingraph[to].emplace_back(from, weight);
// 	}
// 	m_edge_values(from,to) = weight;
	m_processed = false;
}

void DiGraph::add_edge(const std::string& from, const std::string& to, weight_t weight)
{
	assert(m_namemap.find(from) != m_namemap.end());
	assert(m_namemap.find(to) != m_namemap.end());
	add_edge(m_namemap[from], m_namemap[to], weight);
}


void DiGraph::process()
{
	if (!m_processed)
	{
// 		std::cout << "Processing graph..." << std::endl;
		Chronometer C;
		find_weakly_connected_components();
// 		std::cout << "weak coloring = " << m_weak_coloring << std::endl;
// 		std::cout << "before removing, size: " << num_vertices() << std::endl;
		remove_bad_nodes();
// 		std::cout << "After removing, new size: " << m_n << std::endl;

		find_strongly_connected_components();

// 		std::cout << "Scc: " << m_scc_coloring << std::endl;

		m_processed = true;
// 		std::cout << "Finished processing graph after " << C.Peek() << " seconds" << std::endl;
	}
}

Path DiGraph::forward_backward_dfs(node_t start)
{
	Path FP = dfs_search_path_forward(start, Options.dfs_time_woimprovement);

	Path BP = FP;

	int numremoves  = std::min(20, int(BP.size() / 2) - 1);

	for (int i = 0; i < numremoves; ++i)
	{
		BP.pop_front();
	}

	dfs_search_path_backward(BP, Options.dfs_time_woimprovement);

	if (BP.value() > FP.value())
	{
		FP = BP;
	}

	return FP;
}

Path DiGraph::dfs_forward_full() const
{
	Path Best;

	//FORWARD SEARCH

	auto numrestarts = Options.dfs_forward_num_starting_nodes;
	auto t = Options.dfs_time_woimprovement;

	if (numrestarts > num_vertices())
		numrestarts = num_vertices();


	for (int i = 0; i < numrestarts; ++i)
	{
		node_t node = m_basic_topological_ordering[i / 2];
		if ((i%2) == 1)
		{
			node = rand() % num_vertices();
		}

		Path P = dfs_search_path_forward(node, t);

		if (P.value() > Best.value())
		{
			Best = std::move(P);
		}
	}

	int num_deletions = Options.dfs_how_many_to_erase_from_opposite_side;
	auto P = Best;
	int Psize = static_cast<int>(P.size());

	if (2*num_deletions > Psize)
	{
		num_deletions = Psize / 2;
	}

	for (int i = 0; i < num_deletions; ++i)
	{
		P.pop_front();
	}

	dfs_search_path_backward(P, t);

	if (P.value() > Best.value())
	{
		Best = std::move(P);
	}
	
	return Best;
}

Path DiGraph::dfs_backward_full() const
{
	Path Best;

	//FORWARD SEARCH
	
	auto numrestarts = Options.dfs_backward_num_starting_nodes;
	auto t = Options.dfs_time_woimprovement;

	if (numrestarts > num_vertices())
		numrestarts = num_vertices(); //really no point in starting from more than the number of vertices.
	
	for (int i = 0; i < numrestarts; ++i)
	{
		node_t node = m_basic_topological_ordering_in[i/2];
		if ((i%2) == 1)
		{
			node = rand() % num_vertices();
		}

		Path P = dfs_search_path_backward(node, t);

		if (P.value() > Best.value())
		{
			Best = std::move(P);
		}
	}

	int num_deletions = Options.dfs_how_many_to_erase_from_opposite_side;
	auto P = Best;
	int Psize = static_cast<int>(P.size());

	if (num_deletions > Psize/2)
	{
		num_deletions = Psize/2;
	}

	for (int i = 0; i < num_deletions; ++i)
	{
		P.pop_back();
	}

	dfs_search_path_forward(P, t);

	if (P.value() > Best.value())
		Best = std::move(P);
	
	return Best;
}

/////////////////
/// These are just some parameters I found useful. Probably best to do something smarter. 
/// Check trainer.cpp (hpp) to see a way to find these parameters using genetic algorithms.
////////////
param_t GetParams(int i)
{
	if (i == 0)
		return {{1,4,16,64,1,4,16,64}};
	if (i == 1)
		return {{-43,31,11,58,-4,23,43,45}};
	param_t X;
	for (auto& x : X)
		x = rand()%100-10;
	return X;
}

Path DiGraph::dfs_search()
{
	Path Best;
	
	for (int i = 0; i < Options.dfs_num_parameter_restarts; ++i)
	{
		set_parameters(GetParams(i));
		auto ForwardPath = dfs_forward_full();
		
		if (ForwardPath.value() > Best.value())
			Best = std::move(ForwardPath);
		
		auto BackwardPath = dfs_backward_full();
		if (BackwardPath.value() > Best.value())
			Best = std::move(BackwardPath);
	}
	return Best;
}


Path DiGraph::FindLongestSimplePath()
{
	Chronometer C;
	process();

	std::cout << "Doing DFS search..." << std::endl;
	Path best = dfs_search();

	std::cout << "Done DFS search! Value = " << best.value() << std::endl;
	std::cout << "Time taken for dfs: " << C.Peek() << std::endl;

	std::cout << "Doing PTO improving search..." << std::endl;
	pto_search(best);

	return best;
}

void DiGraph::pto_search(Path& A) const
{
	Chronometer C;

	PseudoTopoOrder PTO = get_random_pseudotopological_order();
// 	std::cout << "Before applying A, PTO.Value() = " << PTO.Value() << std::endl;

	PTO.apply(A);

// 	std::cout << "After applying A, PTO.Value() = " << PTO.Value() << std::endl;

	PTO.open_edges_until_no_more_improvement_found();

// 	std::cout << "After opening edges, PTO.Value() = " << PTO.Value() << std::endl;

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
	std::vector<node_t> toRemove;
	toRemove.reserve(m_n / 2);

// 		std::cout << "removing bad nodes!" << std::endl;
// 	int val = dfs_search_path(0.1).Value();

	size_t val = 0;

	for (auto& m_weakly_connected_component : m_weakly_connected_components)
	{
		if (val < m_weakly_connected_component.size())
		{
			val = m_weakly_connected_component.size();
		}
	}

	for (auto& m_weakly_connected_component : m_weakly_connected_components)
	{
		if (m_weakly_connected_component.size() < val)
		{
			for (auto x : m_weakly_connected_component)
			{
				toRemove.push_back(x);
			}
		}
	}

// 	std::cout << "done adding toremove! sorting..." << std::endl;
// 	sort(toRemove.begin(), toRemove.end());

// 	for (int i = 0; i < n; ++i)
// 	{
// 		if (outgraph[i].size() == 0 && ingraph[i].size() == 0)
// 			toRemove.push_back(i);
// 	}
// 	std::cout << "Removing: " << toRemove.size() << " nodes" << std::endl;
	remove_nodes(toRemove);
}


DiGraph DiGraph::with_nodes_removed(std::vector<node_t>& toRemove) const
{
	if (toRemove.empty())
	{
		return *this;
	}

	sort(toRemove.begin(), toRemove.end());

	size_t new_n = m_n - toRemove.size();
// 	DiGraph D(new_n);
	std::vector<std::string> new_names;
	new_names.reserve(new_n);

	std::vector<node_t> removalfunction(m_n, INVALID_NODE);
	std::vector<node_t> removalfunctioninverse(new_n, 0);

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
		}
		else
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
			D.add_edge(v, newneigh, oldneigh.Weight());
		}
	}

	return D;
}

void DiGraph::remove_nodes(std::vector<node_t>& toRemove)
{
	*this = with_nodes_removed(toRemove);
}

void DiGraph::dfs_search_path_forward(Path& P, double maxnumseconds) const
{
	ExpandGreedyBack(*this, P);
// 	std::cout << "(" << ChronometerPeek() << ", " << P.Value() << ")," << std::endl;

	Chronometer C;
	Path Q = P;
	auto comp = [this](node_t a, node_t b)
	{
		return ex_compare(a, b);
	};

	while (C.Peek() < maxnumseconds && dfs_outnext(*this, Q, comp))
	{
		if (Q.value() > P.value())
		{
			P = Q;
			C.Reset();
// 			std::cout << "(" << C.Peek() << ", " << P.value() << ")," << std::endl;
// 			std::cout << P.value()/2 << std::endl;
		}
	}
}

void DiGraph::dfs_search_path_backward(Path& P, double maxnumseconds) const
{
	ExpandGreedyBack(*this, P);
// 	std::cout << "(" << ChronometerPeek() << ", " << P.Value() << ")," << std::endl;

	Path Q = P;
	auto comp = [this](node_t a, node_t b)
	{
		return in_compare(a, b);
	};

	Chronometer C;
	while (C.Peek() < maxnumseconds && dfs_innext(*this, Q, comp))
	{
		if (Q.value() > P.value())
		{
			P = Q;
			C.Reset();
// 			std::cout << "(" << ChronometerPeek() << ", " << P.Value() << ")," << std::endl;
// 			std::cout << P.value()/2 << std::endl;
		}
	}
}

Path DiGraph::dfs_search_path_forward(node_t start, double maxnumseconds) const
{
	Path P(num_vertices(), start);
	dfs_search_path_forward(P, maxnumseconds);
	return P;
}

Path DiGraph::dfs_search_path_backward(node_t start, double maxnumseconds) const
{
	Path P(num_vertices(), start);
	dfs_search_path_backward(P, maxnumseconds);
	return P;
}

void DiGraph::heuristic_processing()
{
	m_heuristic_out.resize(m_n, 0);
	m_heuristic_in.resize(m_n, 0);

	// 	double maxtime = 0.2/n;

// 	#pragma omp parallel for
	for (node_t i = 0; i < m_n; ++i)
	{
		m_heuristic_out[i] = get_heuristic_out(i);
		m_heuristic_in[i] = get_heuristic_in(i);
	}


	m_basic_topological_ordering = range<node_t>(m_n);
	sort(m_basic_topological_ordering.begin(), m_basic_topological_ordering.end(), [this](node_t a, node_t b) -> bool
	{
		//First, order by rank
		if (rank_out(a) < rank_out(b))
		{
			return false;
		}
		if (rank_out(a) > rank_out(b))
		{
			return true;
		}

		//If someone doesn't have out-neighbours, they should go last
		if (outdegree(a) == 0)
		{
			return false;
		}
		if (outdegree(b) == 0)
		{
			return true;
		}

		//The ones who have NO OTHER entry points should go first
		if (indegree(b) == 1)
		{
			return false;
		}
		if (indegree(a) == 1)
		{
			return true;
		}

		//If all else fails, sort by heuristic
		return m_heuristic_out[a] < m_heuristic_out[b];
	});

	m_basic_topological_ordering_inverse.resize(m_n);

	for (node_t i = 0; i < m_n; ++i)
	{
		m_basic_topological_ordering_inverse[m_basic_topological_ordering[i]] = i;
	}

	m_basic_topological_ordering_in = range<node_t>(m_n);
	sort(m_basic_topological_ordering_in.begin(), m_basic_topological_ordering_in.end(), [this](node_t a, node_t b) -> bool
	{
		if (rank_in(a) > rank_in(b))
		{
			return false;
		}
		if (rank_in(a) < rank_in(b))
		{
			return true;
		}


		if (indegree(a) == 0)
		{
			return false;
		}
		if (indegree(b) == 0)
		{
			return true;
		}

		if (outdegree(b) == 1)
		{
			return false;
		}
		if (outdegree(a) == 1)
		{
			return true;
		}

		return m_heuristic_in[a] < m_heuristic_in[b];
	});

	m_basic_topological_ordering_inverse_in.resize(m_n);

	for (node_t i = 0; i < m_n; ++i)
	{
		m_basic_topological_ordering_inverse_in[m_basic_topological_ordering_in[i]] = i;
	}


	for (size_t i = 0; i < m_n; ++i)
	{
		sort(m_outgraph[i].begin(), m_outgraph[i].end(), [this](node_t a, node_t b) -> bool
		{
			return ex_compare(a, b);
		});

		sort(m_ingraph[i].begin(), m_ingraph[i].end(), [this](node_t a, node_t b) -> bool
		{
			return in_compare(a, b);
		});
	}

}

void DiGraph::DFSUtil(node_t v, std::vector<bool>& visited)
{
	visited[v] = true;

	for (auto i : m_outgraph[v])
	{
		if (!visited[i])
		{
			DFSUtil(i, visited);
		}
	}
}


void DiGraph::DFSUtilReversed(node_t v, std::vector< char >& visited, int current)
{
	visited[v] = 1;
// 	std::cout << " v = " << v << " and current = " << current << " and scc.size() = " << strongly_connected_components.size() << std::endl;
	m_strongly_connected_components[current].push_back(v);
	m_scc_coloring[v] = current;

// 	std::cout << "good" << std::endl;
	for (auto i : m_ingraph[v])
	{
		if (visited[i] == 0)
		{
			DFSUtilReversed(i, visited, current);
		}
	}
}


void DiGraph::topo_fill_order(node_t v, std::vector< char >& visited, std::stack< node_t >& Stack)
{
	// Mark the current node as visited and print it
	visited[v] = 1;

	// Recur for all the vertices outgraphacent to this vertex
	for (auto i : m_outgraph[v])
	{
		if (visited[i] == 0)
		{
			topo_fill_order(i, visited, Stack);
		}
	}

	// All vertices reachable from v are processed by now, push v
	Stack.push(v);
}

void DiGraph::find_strongly_connected_components()
{
	m_strongly_connected_components.clear();
	m_strongly_connected_components.reserve(m_n / 4);

	m_scc_coloring.resize(m_n);

	std::stack<node_t> Stack;

	std::vector<char> visited(m_n, 0);

	for (int i = 0; i < m_n; ++i)
	{
		if (!static_cast<bool>(visited[i]))
		{
			topo_fill_order(i, visited, Stack);
		}
	}

	for (int i = 0; i < m_n; i++)
	{
		visited[i] = 0;
	}

	int current = 0;

// 	std::cout << "before starting, scc.size() = " << m_strongly_connected_components.size() << std::endl;
	while (!Stack.empty())
	{
		int v = Stack.top();
		Stack.pop();

		m_strongly_connected_components.emplace_back();

		if (!static_cast<bool>(visited[v]))
		{
			DFSUtilReversed(v, visited, current);
// 			std::cout << std::endl;
			++current;
			m_strongly_connected_components.emplace_back();
		}
	}

//     std::cout << "Finished first loop" << std::endl;

	while (m_strongly_connected_components.back().empty())
	{
		m_strongly_connected_components.pop_back();
	}

//     std::cout << "Finished popping the empty components. WTF?" << std::endl;

	int num_scc = m_strongly_connected_components.size();
	m_scc_rank_in.resize(m_n, 0);

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
				{
					continue;
				}

				int candidate = m_scc_rank_in[connected_component_neigh] + 1;

				if (candidate > m_scc_rank_in[scc])
				{
					m_scc_rank_in[scc] = candidate;
				}
			}
		}
	}

// 	std::cout << "Finished with this weird loop." << std::endl;

	m_scc_rank_out.resize(num_scc, 0);

	for (int scc = num_scc - 1; scc >= 0; --scc)
	{
		for (int x : m_strongly_connected_components[scc])
		{
			// x is the actual node in the scc component
			for (auto v : m_outgraph[x])
			{
				//v is the actual node in the scc component
				int connected_component_neigh = m_scc_coloring[v];

				if (connected_component_neigh == scc)
				{
					continue;
				}

				int candidate = m_scc_rank_out[connected_component_neigh] + 1;

				if (candidate > m_scc_rank_out[scc])
				{
					m_scc_rank_out[scc] = candidate;
				}
			}
		}
	}

// 	std::cout << "Finished with this other weird loop." << std::endl;

	int start = 0;

	for (int color = 0; color < num_scc; ++color)
	{
		int s = m_strongly_connected_components[color].size();

		if (s > 1)
		{
			m_scc_big_components.emplace_back(start, start + s, color);
		}

		start += s;
	}

// 	std::cout << "Finished with this third weird loop." << std::endl;

// 	for (int i = 0; i < m_scc_big_components.size(); ++i)
// 	{
// 		auto u = m_scc_big_components[i];
// 		std::cout << i << ": " << u.start << " --> " << u.end << " w color " << u.color << std::endl;
// 	}

// 	std::cout << "Finished printing big component info." << std::endl;
// 	std::cout << " with m_strongly_connected_components.size() = " << m_strongly_connected_components.size() << std::endl;
// 	std::cout << " and num_scc = " << num_scc << std::endl;
// 	for (int i = 0; i < num_scc; ++i)
// 	{
// 		std::cout << i << std::endl;
// 		auto x = m_strongly_connected_components[i];
// 		if (x.size() > 1)
// 			std::cout << i << ": " << x << std::endl;
// 	}

// 	std::cout << "out of..." << m_strongly_connected_components.size();
// 	std::cout << "scc_rank_out = " << m_scc_rank_out << std::endl;

}

PseudoTopoOrder DiGraph::get_random_pseudotopological_order() const
{
	std::vector<node_t> topo_sort(m_n, 0); // at index i, it has a node topo_sort[i]
	std::vector<node_t> topo_sort_inverse(m_n, 0); // at node i, x=topo_sort_inverse[i] is such that topo_sort[x] = i
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

	return PseudoTopoOrder(*this, topo_sort, topo_sort_inverse);
}

void DiGraph::DFSUtilWeak(node_t node, int color)
{
	m_weak_coloring[node] = color;
	m_weakly_connected_components.back().push_back(node);

// 	m_weakly_connected_components_values.back() += vertex_values[node]-1;
	for (auto x : m_outgraph[node])
	{
		if (m_weak_coloring[x] != -1)
		{
			continue;
		}

		DFSUtilWeak(x, color);
	}

	for (auto x : m_ingraph[node])
	{
		if (m_weak_coloring[x] != -1)
		{
			continue;
		}

		DFSUtilWeak(x, color);
	}
}

void DiGraph::find_weakly_connected_components()
{
	m_weak_coloring.clear();
	m_weak_coloring.resize(m_n, -1);
	m_weakly_connected_components.clear();
// 	m_weakly_connected_components_values.clear();
	int minvalidcoloring = 0;

	for (int start = 0; start < m_n; ++start)
	{
		if (m_weak_coloring[start] != -1)
		{
			continue;
		}

		m_weakly_connected_components.emplace_back();
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
				heuristicex += a3 + a4 * m_outgraph[z].size();
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
			heuristicin  += a2 + a3 * m_ingraph[y].size();

			for (auto z : m_ingraph[y])
			{
				heuristicin  += a3 + a4 * m_ingraph[z].size();
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
		for (int j = i + 1; j < n; ++j)
		{
			if (probability_of_true(p))
			{
				int a = i;
				int b = j;

				if ((rand()%2) == 1)
				{
					std::swap(a, b);
				}

				D.add_edge(a, b);
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
		for (int j = i + 1; j < n; ++j)
		{
			if (probability_of_true(p))
			{
				int a = i;
				int b = j;
				weight_t w = random_real(minweight, maxweight);

				if (rand() % 2 == 1)
				{
					std::swap(a, b);
				}

				D.add_edge(a, b, w);
			}
		}
	}

// 	D.process();
	return D;
}

std::ostream& operator<<(std::ostream& os, const param_t& a)
{
	int i = 0;
	for (auto x : a)
	{
		os << x << ' ';
		if (i == 4)
			os << '|';
		++i;
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, const DiGraph& M)
{
	os << "Digraph on " << M.num_vertices() << " vertices: " << M.get_vertex_names();

	for (int i = 0; i < M.num_vertices(); ++i)
	{
		const std::string& name = M.get_vertex_name(i);

		for (auto v : M.outneighbors(i))
		{
			os << std::endl;
			const std::string& vname = M.get_vertex_name(v);
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
		{
			break;
		}

		P.emplace_back(t, t.Weight());
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
		{
			break;
		}

		P.emplace_front(t, t.Weight());
	}
}
