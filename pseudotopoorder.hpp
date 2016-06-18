#pragma once

#include "moviegraph.hpp"
#include "path.hpp"

class PseudoTopoOrder
{
public:
	PseudoTopoOrder(const DiGraph& m, const vector<node_t>& ts, const vector<node_t>& tsi) : 
	pto(ts), 
	pto_inverse(tsi), 
	dynamic_programming(ts.size(),0), 
	best_index(0), 
	first_unknown(0),
	path_filled(false),
	m_path(),
	m_parent(m) {}

	size_t size() const { return dynamic_programming.size(); }
	
	inline void set(int i, int v)
	{
		pto[i] = v;
		pto_inverse[v] = i;
		AnnounceModification(i);
	}

	inline void transpose(int i, int j)
	{
		swap(pto[i],pto[j]);
		pto_inverse[pto[i]] = i;
		pto_inverse[pto[j]] = j;
		AnnounceModification(i);
		AnnounceModification(j);
	}
	
	void transfer(int a, int b, int c, int d, int h);
	
	void reverse_order(int a, int b);
	
	Path get_path();

	void shuffle(int  a, int b); // this assumes already a < b and they belong to the same component
	sumweight_t  Value() { FillDP(); return dynamic_programming[best_index]; }
	void apply(const Path& P);
	void apply(const Path& P, int u, int v);
// 	void random_apply(const Path& P);
	
	bool eXtreme_edge_opener();
	
	void open_edges_until_no_more_improvement_found(double maxnumseconds);
	
	void open_edge();
	
	void randomize();
	
	// sorts range from a to b (must be same scc) according to a heuristic so as to maximize the improvement chance.
	void heuristic_sort(int a, int b, int numtimes);
	int get_outneighbor_in_range(int a, int b, node_t node);

private:
	void FillDP();
	void FillPath();
	void RecalcTopoInverse();
	
	inline void AnnounceModification(node_t i) { if  (i < first_unknown) first_unknown = i; path_filled = false; }
	inline void transpose_na(int i, int j)
	{
		swap(pto[i],pto[j]);
		pto_inverse[pto[i]] = i;
		pto_inverse[pto[j]] = j;
	}
	
	inline void set_na(int i, int v)
	{
		pto[i] = v;
		pto_inverse[v] = i;
	}
	
	vector<node_t> pto;
	vector<node_t> pto_inverse;
	vector<sumweight_t> dynamic_programming;
	int best_index;
	int first_unknown;
	bool path_filled;
	vector<NeighborNode> m_path; //Path is filled with indices in REVERSE ORDER
	const DiGraph& m_parent;
};
