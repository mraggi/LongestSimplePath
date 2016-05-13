#include "path.hpp"
#include <cassert>

void Path::AddNodeBack(node_t t)
{
    m_explored[t] = true;
	node_t l = m_path.back();
	m_totalvalue += m_parent->edge_value(l,t);
    m_path.push_back(t);
}

void Path::AddNodeFront(node_t t)
{
    m_explored[t] = true;
    node_t l = m_path.front();
    m_totalvalue += m_parent->edge_value(t,l);
    m_path.push_front(t);
}

void Path::PopBack()
{
// 	assert (m_path.size() > 1 && "Cannot have empty paths!");
	if (m_path.size() < 2) return;
    node_t l = m_path.back();
    m_explored[l] = false;
    m_path.pop_back();
    node_t ntl = m_path.back();
    m_totalvalue -= m_parent->edge_value(ntl,l);
}

void Path::PopFront()
{
	if (m_path.size() < 2) return;
    node_t first = m_path.front();
    m_explored[first] = false;
    m_path.pop_front();
    node_t second = m_path.front();
    m_totalvalue -= m_parent->edge_value(first,second);
}

void Path::ExpandGreedyBack()
{
	while (true)
	{
		auto l = m_path.back();
		auto& Neighs = m_parent->exneighbors(l);
		
		auto t = FirstNotExplored(Neighs);
		if (t == INVALID_NODE)
			break;
		AddNodeBack(t);
	}
}

void Path::transform_nodes(const vector< node_t >& m_removalfunctioninverse)
{
	for (auto& x : m_path)
		x = m_removalfunctioninverse[x];
}


void Path::ExpandGreedyFront()
{
	while (true)
	{
		auto l = m_path.front();
		auto& Neighs = m_parent->inneighbors(l);
		
		auto t = FirstNotExplored(Neighs);
		if (t == INVALID_NODE)
			break;
		AddNodeFront(t);
	}
}

int Path::FirstNotExploredEX(const vector<node_t>& Nodes, node_t start) const
{
	auto it = std::upper_bound(Nodes.begin(), Nodes.end(), start, [this] (node_t a, node_t b) -> bool
		{
			return m_parent->ex_compare(a,b);
		});
// 	++it;
	while (it != Nodes.end() && m_explored[*it])
		++it;
	if (it == Nodes.end())
		return INVALID_NODE;
	return *it;
}


int Path::FirstNotExploredIN(const vector<node_t>& Nodes, node_t start) const
{
	auto it = std::upper_bound(Nodes.begin(), Nodes.end(), start, [this] (node_t a, node_t b) -> bool
		{
			return m_parent->in_compare(a,b);
		});
// 	++it;
	while (it != Nodes.end() && m_explored[*it])
		++it;
	if (it == Nodes.end())
		return INVALID_NODE;
	return *it;
	
}

node_t Path::FirstNotExplored(const vector<node_t>& Nodes, node_t start) const
{
	bool seenstart = false;
	for (auto x : Nodes)
	{
		if (x == start)
		{
			seenstart = true;
			continue;
		}	
		
		if (seenstart && m_explored[x] == false)
			return x;
	}
	return INVALID_NODE;
}

bool Path::DFSBackNext(int maxbacktrack)
{
	int lastNode = m_path.back();
    const vector<node_t>* Neighs = &(m_parent->exneighbors(lastNode));

	int t = FirstNotExplored(*Neighs);

    while (t == -1 && m_path.size() > 1 && maxbacktrack > 0) //this means all nodes in Neigh have been explored
    {
		--maxbacktrack;
		lastNode = m_path.back();
		PopBack();
		int father = m_path.back();
		Neighs = &(m_parent->exneighbors(father));
		t = FirstNotExploredEX(*Neighs,lastNode);
    }
    if (t == -1)
		return false; // this means we have finished DFS!!
	AddNodeBack(t);
	ExpandGreedyBack();
    return true;
}


bool Path::DFSFrontNext(int maxbacktrack)
{
	int lastNode = m_path.front();
// 	cout << "lastnode = " << lastNode << endl;
    const vector<node_t>* Neighs = &(m_parent->inneighbors(lastNode));
// 	cout << "Neighs = " << *Neighs << endl;

	int t = FirstNotExplored(*Neighs);
// 	cout << "t = " << t << endl;

    while (t == -1 && m_path.size() > 1 && maxbacktrack > 0) //this means all nodes in Neigh have been explored
    {
		--maxbacktrack;
		lastNode = m_path.front();
		PopFront();
		int father = m_path.front();
		Neighs = &(m_parent->inneighbors(father));
		t = FirstNotExploredIN(*Neighs,lastNode);
    }
    if (t == -1)
		return false;
	AddNodeFront(t);
	ExpandGreedyFront();
    return true;
}

Path::Path(const DiGraph* const parent, const deque<node_t>& path) : 	m_parent(parent),
											m_totalvalue(0),
											m_path(path),
											m_explored(m_parent->get_size(),0)
{
	node_t prev = m_path.front();
	for (auto x : m_path)
	{
		m_explored[x] = true;
		m_totalvalue += m_parent->edge_value(prev,x);
		prev = x;
	}
}

Path::Path(const DiGraph* const parent, const deque<node_t>& path, sumweight_t value) : 	m_parent(parent),
																					m_totalvalue(value),
																					m_path(path),
																					m_explored(m_parent->get_size(),0)
{
	for (auto x : m_path)
	{
		m_explored[x] = true;
	}
}

void Path::Reset()
{
	for (auto x : m_path)
	{
		m_explored[x] = 0;
	}
	m_path.clear();
	m_totalvalue = 0;

}


Path Path::BestCompletion()
{
	Path best = *this;
	vector<Path> frontier = {*this};
	while (!frontier.empty())
	{
		Path P = frontier.back();
		frontier.pop_back();
		if (P.m_totalvalue > best.m_totalvalue)
			best = P;
		
		int lastnode = P.m_path.back();
		auto exx = m_parent->exneighbors(lastnode);
		random_shuffle(exx.begin(), exx.end());
		
		for (auto x : exx)
		{
			if (P.m_explored[x] == false)
			{
				
				Path PP = P;
				PP.AddNodeBack(x);
				frontier.push_back(PP);
				
			}
		}
	}
	return best;
}

bool Path::operator==(const Path& P) const
{
	return (m_totalvalue == P.m_totalvalue && m_path == P.m_path);
}

bool Path::operator!=(const Path& P) const
{
	return !(*this == P);
}

std::ostream& operator<<(std::ostream& os, const Path& P)
{
    auto B = P.get_path();
    for (size_t i = 0; i < B.size()-1; ++i)
	{
		os << P.get_digraph()->get_vertex_name(B[i]) << u8" ⟼ ";
	}
	os << P.get_digraph()->get_vertex_name(B.back()) << endl;
	return os;
}
