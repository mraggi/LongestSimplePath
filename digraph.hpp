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

class DiGraph;
class Path;
class PseudoTopoOrder;

using namespace std;

struct SccInfo
{
	SccInfo(int s, int e, int c) : start(s), end(e), color(c) {}
	int start;
	int end;
	int color;
};

typedef int node_t;
typedef char weight_t; // Do NOT change this to float unless you know what you are doing: little errors give the impression it's not working.
typedef int sumweight_t;

using ParamType = std::array<double, 8>;

const node_t INVALID_NODE = -1; //make sure it's supported and not a "common" node (like 0, 1, 2...)

class DiGraph
{
public:
    explicit DiGraph(node_t numberOfNodes);
    explicit DiGraph(const vector<string>& vertex_names);

	//		Graph modification functions
	void add_edge(node_t from, node_t to, weight_t weight = 1);
	void add_edge(const string& from, const string& to, weight_t weight = 1);
	
	// Find connected components, heuristics, etc.
	void process(); // WARNING! Since it removes and renames nodes, after "processing" your nodes might have been renamed!

	// Get Graph Info
	node_t get_size() const { return m_n; }
	node_t num_vertices() const { return m_n; }
	size_t num_edges() const;
    
    const string& get_vertex_name(node_t i) const { return m_node_names[i]; }
    node_t get_vertex_index(const string& name) const 
    { 
        auto it = m_namemap.find(name);
        return (*it).second; 
        
    }
    const vector<string>& get_vertex_names() const { return m_node_names; } 
    
    void set_parameters(const ParamType& new_params) 
    {
        m_params = new_params;
        heuristic_processing();
    }
    
	int rank_ex(node_t node) const { return m_scc_rank_ex[m_scc_coloring[node]]; }
	int rank_in(node_t node) const { return m_scc_rank_in[m_scc_coloring[node]]; }
	
	inline const vector<node_t>& exneighbors(node_t n) const { return m_exgraph[n]; }
	inline const vector<node_t>& inneighbors(node_t n) const { return m_ingraph[n]; }
	inline const weight_t edge_value(node_t from, node_t to) const { return m_edge_values(from,to); }
	inline const vector<vector<node_t>> strongly_connected_components() const { return m_strongly_connected_components; }
	inline const vector<SccInfo> big_scc() const { return m_scc_big_components; }
	
	// This is the order in which the exneighbors are sorted
	inline bool ex_compare(node_t a, node_t b) const { return m_basic_topological_ordering_inverse[a] < m_basic_topological_ordering_inverse[b]; }

	// This is the order in which the exneighbors are sorted
	inline bool in_compare(node_t a, node_t b) const { return m_basic_topological_ordering_inverse_in[a] < m_basic_topological_ordering_inverse_in[b]; }
	
	// Functions related to paths
	
	void dfs_search_path_forward(Path& P, double maxnumseconds) const;
	void dfs_search_path_reverse(Path& P, double maxnumseconds) const;
	
	Path dfs_search_path_forward(node_t start, double maxnumseconds) const;
	Path dfs_search_path_reverse(node_t start, double maxnumseconds) const;

	Path dfs_search(double maxnumsecondswithoutimprovement) const;
    void pto_search(Path& A, double maxnumseconds) const;
    
	// gets a random path, by shuffling the exgraph and ingraph orders
	Path get_random_path(double maxnumseconds) const;
	
	PseudoTopoOrder get_random_pseudotopological_order() const;
	//Paths
	Path FindLongestSimplePath(double numseconds);
	Path FindLongestSimplePathPureDFS(double numseconds);
	
	bool TopologicalLessThan(node_t a, node_t b) const { return m_basic_topological_ordering_inverse[a] < m_basic_topological_ordering_inverse[b]; }
	
	static DiGraph CreateRandomDiGraph(int n, double p);
	static DiGraph CreateRandomWeightedDiGraph(int n, double p, weight_t minweight, weight_t maxweight);
	
private:
    // Utils for creating the graph
	void remove_bad_nodes();
	void remove_nodes(const vector<node_t>& toRemove);
	void heuristic_processing();
	double get_heuristic_ex(node_t node);
	double get_heuristic_in(node_t node);
	
// 	void branch_and_bound();
	
	// Utils to find connected components
	void find_weakly_connected_components();
	void find_strongly_connected_components();
	void topo_fill_order( node_t v, vector< char >& visited, stack< node_t >& Stack ); 
    void DFSUtil( node_t v, vector< bool >& visited );
    void DFSUtilReversed( node_t v, vector< char >& visited, int current );
    void DFSUtilWeak(node_t start, int minvalidcoloring);
	
protected:
	// DiGraph insides
	node_t m_n;
	vector<vector<node_t>> m_exgraph;
	vector<vector<node_t>> m_ingraph;
	SquareSparseMatrix<weight_t> m_edge_values;

private:
	bool m_processed;
	// Connected Components
	vector<vector<node_t>> m_strongly_connected_components;
	vector<node_t> m_scc_coloring;
	vector<node_t> m_scc_rank_ex; // This is the ex rank of the connected components
	vector<node_t> m_scc_rank_in; // This is the in rank of the connected components
	
	vector<short> m_weak_coloring;
	vector<vector<node_t>> m_weakly_connected_components;

	//Heuristics
	vector<double> m_heuristic_ex;
	vector<double> m_heuristic_in;
	
	vector<node_t> m_basic_topological_ordering;
	vector<node_t> m_basic_topological_ordering_in;
	vector<node_t> m_basic_topological_ordering_inverse;
	vector<node_t> m_basic_topological_ordering_inverse_in;
	
	vector<SccInfo> m_scc_big_components;
	
    vector<string> m_node_names;
    unordered_map<string, node_t> m_namemap;
    
    
    ParamType m_params {{-43,31,11,58,-4,23,43,45}};
    
	friend class PseudoTopoOrder;
};

std::ostream& operator<<(std::ostream& os, const DiGraph& M);
std::ostream& operator<<(std::ostream& os, const ParamType& a);
