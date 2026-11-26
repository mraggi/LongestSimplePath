#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <ctime>
#include <stack>
#include <algorithm>
#include <set>
#include "helpers.hpp"
#include "path.hpp"


struct SccInfo
{
	SccInfo(int s, int e, int c) : start(s), end(e), color(c) {}
	int start;
	int end;
	int color;
};
class PseudoTopoOrder;

using param_t = std::array<double, 8>;

class DiGraph
{
public:
	explicit DiGraph(node_t numNodes);
	
	DiGraph(const std::vector<std::string>& vertex_names);

	//		Graph modification functions
	void add_edge(node_t from, node_t to, weight_t weight = 1);
	void add_edge(const std::string& from, const std::string& to, weight_t weight = 1);

	// Find connected components, heuristics, etc.
	void process(); // WARNING! Since it removes and renames nodes, after "processing" your nodes might have been renamed!

	// Get Graph Info
	node_t get_size() const { return m_n; }
	node_t num_vertices() const { return m_n; }
	size_t num_edges() const;
	size_t outdegree(node_t node) const { return m_outgraph[node].size(); }
	size_t indegree(node_t node) const { return m_ingraph[node].size(); }

	bool are_neighbors(node_t a, node_t b) const;
	
	const std::string& get_vertex_name(node_t i) const { return m_node_names[i]; }
	node_t get_vertex_index(const std::string& name) const
	{
		auto it = m_namemap.find(name);
		return (*it).second;

	}
	const std::vector<std::string>& get_vertex_names() const { return m_node_names; }

	bool check_path(const Path& P) const;
	
	void set_parameters(const param_t& new_params)
	{
		m_params = new_params;
		heuristic_processing();
	}

	int rank_out(node_t node) const { return m_scc_rank_out[m_scc_coloring[node]]; }
	int rank_in(node_t node) const { return m_scc_rank_in[m_scc_coloring[node]]; }

	inline const std::vector<NeighborNode>& outneighbors(node_t n) const { return m_outgraph[n]; }
	inline const std::vector<NeighborNode>& inneighbors(node_t n) const { return m_ingraph[n]; }
// 	inline const weight_t edge_value(node_t from, node_t to) const { return m_edge_values(from,to); }
	inline const std::vector<std::vector<node_t>> strongly_connected_components() const { return m_strongly_connected_components; }
	inline const std::vector<SccInfo> big_scc() const { return m_scc_big_components; }

	// This is the order in which the outneighbors are sorted
	inline bool ex_compare(node_t a, node_t b) const { return m_basic_topological_ordering_inverse[a] < m_basic_topological_ordering_inverse[b]; }

	// This is the order in which the outneighbors are sorted
	inline bool in_compare(node_t a, node_t b) const { return m_basic_topological_ordering_inverse_in[a] < m_basic_topological_ordering_inverse_in[b]; }

	// Functions related to paths

	void dfs_search_path_forward(Path& P, double maxnumseconds) const;
	void dfs_search_path_backward(Path& P, double maxnumseconds) const;

	Path dfs_search_path_forward(node_t start, double maxnumseconds) const;
	Path dfs_search_path_backward(node_t start, double maxnumseconds) const;

	Path dfs_forward_full() const;
	Path dfs_backward_full() const;
	
	Path dfs_search();
	void pto_search(Path& A) const;

	PseudoTopoOrder get_random_pseudotopological_order() const;


	struct SearchOptions
	{
		// DFS part
		int dfs_num_parameter_restarts{3};
		double dfs_time_woimprovement{0.05};
		int dfs_forward_num_starting_nodes{3};
		int dfs_backward_num_starting_nodes{1};
		
		int dfs_how_many_to_erase_from_opposite_side{20};
// 		int num_saved_for_pto {1};

		// PSO part
		double pto_time_without_improvement{2.0};
		int pto_num_times_restart {1};
		int pto_num_heuristic_sort{10000};
		int pto_scc_size_max_pointless {4};
	};
	
	SearchOptions Options {};

	//Paths
	Path FindLongestSimplePath();

	Path FindLongestSimplePathPureDFS();

	bool TopologicalLessThan(node_t a, node_t b) const { return m_basic_topological_ordering_inverse[a] < m_basic_topological_ordering_inverse[b]; }

	static DiGraph CreateRandomDiGraph(int n, double p);
	static DiGraph CreateRandomWeightedDiGraph(int n, double p, weight_t minweight, weight_t maxweight);

	std::vector<DiGraph> StronglyConnectedInducedGraphs();
	
private:
	// Utils for creating the graph
	void remove_bad_nodes();
	void remove_nodes(std::vector<node_t>& toRemove);

	DiGraph with_nodes_removed(std::vector<node_t>& toRemove) const;

	void heuristic_processing();
	double get_heuristic_out(node_t node);
	double get_heuristic_in(node_t node);

// 	void branch_and_bound();

	// Utils to find connected components
	void find_weakly_connected_components();
	void find_strongly_connected_components();
	void topo_fill_order(node_t v, std::vector< char >& visited, std::stack< node_t >& Stack);
	void DFSUtil(node_t v, std::vector< bool >& visited);
	void DFSUtilReversed(node_t v, std::vector< char >& visited, int current);
	void DFSUtilWeak(node_t node, int color);

protected:
	// DiGraph insides
	node_t m_n{};
	std::vector<std::vector<NeighborNode>> m_outgraph;
	std::vector<std::vector<NeighborNode>> m_ingraph;

private:
	Path forward_backward_dfs(node_t start);

	bool m_processed{};
	// Connected Components
	std::vector<std::vector<node_t>> m_strongly_connected_components;
	std::vector<node_t> m_scc_coloring;
	std::vector<node_t> m_scc_rank_out; // This is the ex rank of the connected components
	std::vector<node_t> m_scc_rank_in; // This is the in rank of the connected components

	std::vector<short> m_weak_coloring;
	std::vector<std::vector<node_t>> m_weakly_connected_components;

	//Heuristics
	std::vector<double> m_heuristic_out;
	std::vector<double> m_heuristic_in;

	std::vector<node_t> m_basic_topological_ordering;
	std::vector<node_t> m_basic_topological_ordering_in;
	std::vector<node_t> m_basic_topological_ordering_inverse;
	std::vector<node_t> m_basic_topological_ordering_inverse_in;

	std::vector<SccInfo> m_scc_big_components;

	std::vector<std::string> m_node_names;
	std::unordered_map<std::string, node_t> m_namemap;


//     ParamType m_params {{-43,31,11,58,-4,23,43,45}};
	param_t m_params {{1, 4, 16, 64, 1, 4, 16, 64}};

	friend class PseudoTopoOrder;
	
	mutable sumweight_t global_best {-1};
	mutable Chronometer global_chrono {};
};

std::ostream& operator<<(std::ostream& os, const DiGraph& M);
std::ostream& operator<<(std::ostream& os, const param_t& a);



void ExpandGreedyBack(const DiGraph& G, Path& P);


void ExpandGreedyFront(const DiGraph& G, Path& P);

template <class Compare>
bool dfs_outnext(const DiGraph& G, Path& P, Compare comp)
{
	auto lastNode = P.back();

	auto Neighs = &G.outneighbors(lastNode);

	auto t = P.first_not_explored(*Neighs);

	while (t == INVALID_NODE && P.size() > 1) //this means all nodes in Neigh have been explored
	{
		lastNode = P.back();
		P.pop_back();
		int father = P.back();
		Neighs = &G.outneighbors(father);
		t = P.first_not_explored_binary(*Neighs, lastNode, comp);
	}

	if (t == INVALID_NODE)
		return false; // this means we have finished DFS!!

	P.push_back(t);
	ExpandGreedyBack(G, P);
	return true;
}

template <class Compare>
bool dfs_innext(const DiGraph& G, Path& P, Compare comp)
{
	auto firstNode = P.front();

	auto Neighs = &G.inneighbors(firstNode);

	auto t = P.first_not_explored(*Neighs);

	while (t == INVALID_NODE && P.size() > 1) //this means all nodes in Neigh have been explored
	{
		firstNode = P.front();
		P.pop_front();
		int father = P.front();
		Neighs = &G.inneighbors(father);
		t = P.first_not_explored_binary(*Neighs, firstNode, comp);
	}

	if (t == INVALID_NODE)
		return false; // this means we have finished DFS!!

	P.push_front(t);
	ExpandGreedyFront(G, P);
	return true;
}

// inline bool dfs_outnext(const DiGraph& G, Path& P)
// {
//     return dfs_outnext(G,P,std::less<node_t>());
// }
//
// inline bool dfs_innext(const DiGraph& G, Path& P)
// {
//     return dfs_innext(G,P,std::less<node_t>());
// }
