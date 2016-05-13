#pragma once
#include "digraph.hpp"
#include <deque>

class Path
{
public:
    Path(const DiGraph* parent, node_t node) :m_parent(parent),
													m_totalvalue(0),
													m_path(1,node),
													m_explored(parent->get_size(),0)
	{
		m_explored[node] = 1;
	}
	
	Path(const DiGraph*const parent, const deque<node_t>& path);
	Path(const DiGraph*const parent, const deque<node_t>& path, sumweight_t value);
	
    const DiGraph* const get_digraph() const { return m_parent; }
    
	bool DFSBackNext(int maxbacktrack = 99999);
	bool DFSFrontNext(int maxbacktrack = 99999);
	
	void AddNodeBack(node_t t);
	void AddNodeFront(node_t t);
	
	void PopBack();
	void PopFront();
	void printMovies();
	void Reset();
	
	sumweight_t Value() const { return m_totalvalue; }
	
	void ExpandGreedyBack();
	void ExpandGreedyFront();
	
	node_t FirstNotExplored(const vector<node_t>& Nodes) const
	{
		for (auto x : Nodes)
		{
			if (m_explored[x] == false)
				return x;
		}
		return INVALID_NODE;
	}
	
	bool SanityCheck() const;
	
	node_t FirstNotExplored(const vector< node_t >& Nodes, node_t start) const;
	int FirstNotExploredEX(const vector< node_t >& Nodes, node_t start) const;
	int FirstNotExploredIN(const vector< node_t >& Nodes, node_t start) const;
// 	int LengthOfFirstRunNI(const Path& P) const;
	
	// destroys path
	Path BestCompletion();
	
	bool operator==(const Path& P) const;
	bool operator!=(const Path& P) const;
	const deque<node_t>& get_path() const { return m_path; }
	bool IsNodeInPath(node_t node) const { return m_explored[node]; }
    void transform_nodes(const vector< node_t >& m_removalfunctioninverse);
	
private:
	const DiGraph* m_parent;
    sumweight_t m_totalvalue;
	deque<node_t> m_path;
	vector<char> m_explored;
};

std::ostream& operator<<(std::ostream& os, const Path& P);

