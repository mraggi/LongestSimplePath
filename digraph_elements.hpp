#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <array>
#include <string>
#include <queue>
#include <algorithm>
#include <cassert>
using namespace std;

using node_t = int;
const node_t INVALID_NODE = -1;

using weight_t = int;

// something larger than weight_t for when you have that weight_t doesn't properly hold a sum of weight_t
using sumweight_t = int;
const sumweight_t INF = 200000000;

struct NeighborNode
{
	explicit NeighborNode() : node(INVALID_NODE), weight(0) {} 
	
	explicit NeighborNode(node_t v, weight_t w = 1) : node(v), weight(w) {}
	
	inline operator node_t() const
	{
		return node;
	}
	
	weight_t Weight() const
	{
		return weight;
// 		return 1;
	}
	
	node_t node;
	weight_t weight{1}; //comment 
};

class Path
{
public:
    Path(size_t n) : m_path(), m_explored(n), m_value(0) {}
    Path(size_t n, node_t initialnode) : m_path(), m_explored(n), m_value(0) 
    {
        emplace_back(initialnode,0);
    }
    inline operator const deque<NeighborNode>&() const { return m_path; }
//     inline operator vector<NeighborNode>&()  { return m_path; }
    long value() const { return m_value; }
    long cost() const { return m_value; }
    long weight() const { return m_value; }
    
    void push_back(const NeighborNode& v)
    {
        assert(v.weight == 0 || !m_path.empty());
        assert(v < m_explored.size());

        ++m_explored[v];
        m_value += v.weight;
        m_path.push_back(v);
    }
    
    void emplace_back(node_t node, weight_t weight = 1)
    {
        assert(weight == 0 || !m_path.empty());
        assert(node < m_explored.size());
        ++m_explored[node];
        m_value += weight;
        m_path.emplace_back(node,weight);
    }
    
    void push_front(NeighborNode v)
	{
//         swap(m_path.front().weight,v.weight);
        assert(!m_path.empty() && "Use push_back when it's empty");
        assert(v < m_explored.size());
        m_path.front().weight = v.weight;
        m_path.emplace_front(v,0);
		m_value += v.weight;
        ++m_explored[v];
	}
	
    void emplace_front(node_t node, weight_t weight = 1)
	{
//         auto w = m_path.front().weight;
        assert(!m_path.empty() && "Use emplace_back when it's empty");
        assert(node < m_explored.size());
		m_path.front().weight = weight;
		m_path.emplace_front(node,0);
		m_value += weight;
        ++m_explored[node];
	}
    
    void pop_back()
    {
        assert(!m_path.empty() && "Can't pop when it's already empty!");
        auto v = m_path.back();
        --m_explored[v];
        m_value -= v.weight;
        m_path.pop_back();
    }
    
    void pop_front()
    {
        assert(!m_path.empty() && "Can't pop when it's already empty!");
        auto v = m_path.front();
        --m_explored[v];
        m_path.pop_front();
        m_value -= m_path.front().weight;
        m_path.front().weight = 0;
    }
    
    void clear()
    {
        m_value = 0;
//         m_explored = vector<char>(m_explored.size(),0);
        auto n = m_explored.size();
        m_explored.clear();
        m_explored.resize(n,0);
        m_path.clear();
    }
    
    NeighborNode operator[](size_t i) const { return m_path[i]; }
    NeighborNode& operator[](size_t i) { return m_path[i]; }
    
    bool empty() const
    {
        return m_path.empty();
    }
    
    size_t size() const
    {
        return m_path.size();
    }
    
    template <class Compare>
    NeighborNode first_not_explored_binary(const vector<NeighborNode>& Nodes, node_t start, Compare comp) const
    {
        auto it = std::upper_bound(Nodes.begin(), Nodes.end(), start, comp);
    // 	++it;
        while (it != Nodes.end() && m_explored[*it])
            ++it;
        if (it == Nodes.end())
            return NeighborNode(INVALID_NODE);
        return *it;
    }
    
    NeighborNode first_not_explored_binary(const vector<NeighborNode>& Nodes, node_t start) const
    {
        return first_not_explored_binary(Nodes,start,std::less<node_t>());
    }
    
    NeighborNode first_not_explored(const vector<NeighborNode>& Nodes, node_t start) const
    {
        bool seenstart = false;
        for (auto x : Nodes)
        {
            if (x == start)
            {
                seenstart = true;
                continue;
            }	
            
            if (seenstart && !m_explored[x])
                return x;
        }
        return NeighborNode(INVALID_NODE);
    }
    
    NeighborNode first_not_explored(const vector<NeighborNode>& Nodes) const
    {
        for (auto x : Nodes)
        {
            if (!m_explored[x])
                return x;
        }
        return NeighborNode(INVALID_NODE);
    }
    
    NeighborNode back() const
    {
        return m_path.back();
    }
    
    NeighborNode front() const
    {
        return m_path.front();
    }
    
    
private:
    deque<NeighborNode> m_path;
    vector<char> m_explored;
//     vector<bool> m_explored;
    long m_value;
public:
    const decltype(m_path)& data() const { return m_path; }

};

std::ostream& operator<<(std::ostream& os, const Path& P);

struct Edge
{
	Edge() : from(INVALID_NODE), to(INVALID_NODE), weight(0) {}
	Edge(node_t f, node_t t, weight_t w = 1) : from(f), to(t), weight(w) {}
	node_t operator[](bool i)
	{
		if (i)
			return to;
		return from;
	}
	node_t from;
	node_t to;
	weight_t weight;
};
