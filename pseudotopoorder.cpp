#include "pseudotopoorder.hpp"


void PseudoTopoOrder::shuffle(int a, int b)
{
	random_shuffle(pto.begin()+a, pto.begin()+b);
	for (int i = a; i < b; ++i)
	{
		pto_inverse[pto[i]] = i;
	}
	AnnounceModification(a);
}

void PseudoTopoOrder::RecalcTopoInverse()
{
	int n = pto.size();
	for (int i = 0; i < n; ++i)
	{
		pto_inverse[pto[i]] = i;
	}
	AnnounceModification(0);
}

void PseudoTopoOrder::apply(const Path& P)
{
	auto it = P.get_path().begin();
	int fu = -1;
	int n = pto.size();
	for (int i = 0; i < n; ++i)
	{
		int x = pto[i];
		if (P.IsNodeInPath(x))
		{
			if (fu == -1)
				fu = i;
			pto[i] = *it;
			pto_inverse[*it] = i;
			++it;
		}
	}
	AnnounceModification(fu);
}

void PseudoTopoOrder::apply(const Path& P, int u, int v)
{
// 	cout << "Before applying Q, we get: " << Value() << endl;

	vector<int> indexesofPbetweenuandv;
	vector<int> nodesofPbetweenuandv;
	for (auto x : P.get_path())
	{
		int index = pto_inverse[x];
		if (u <= index && index < v)
		{
			indexesofPbetweenuandv.push_back(index);
			nodesofPbetweenuandv.push_back(x);
// 			cout << "found " << x << " at " << index << endl;
		}
	}
	sort(indexesofPbetweenuandv.begin(), indexesofPbetweenuandv.end());
	int i = 0;
	for (auto p : indexesofPbetweenuandv)
	{
		pto[p] = nodesofPbetweenuandv[i];
		pto_inverse[pto[p]] = p;
		++i;
	}
// 	AnnounceModification(indexesofPbetweenuandv[0]);
	AnnounceModification(u);
// 	cout << "After applying Q, we get: " << Value() << endl;
}


Path PseudoTopoOrder::get_path()
{
	FillPath();
// 	cout << "Filling path, with value: " << Value() << endl;
	deque<node_t> P;
	for (auto i : m_path)
	{
		P.push_front(pto[i]);
	}
	return Path(&m_parent, P, Value());
}

void PseudoTopoOrder::FillDP()
{
	int n = pto.size();
	int best_val = 0;
	if (best_index < first_unknown)
		best_val = dynamic_programming[best_index];
	
	for ( ; first_unknown < n; ++first_unknown)
	{
		int u = pto[first_unknown];
		dynamic_programming[first_unknown] = 0;
		
		auto& neigh = m_parent.inneighbors(u);
		for (auto v : neigh)
		{
			auto j = pto_inverse[v];
			if (first_unknown < j) // should be ignored
				continue;
			auto candidate = dynamic_programming[j] + m_parent.edge_value(v,u);
			if (candidate > dynamic_programming[first_unknown])
			{
				dynamic_programming[first_unknown] = candidate;
				if (candidate > best_val)
				{
// 					cout << "Fill: new best value at " << i << " with value " << best_val << endl;
					best_val = candidate;
					best_index = first_unknown;
				}
			}
		}
	}
// 	cout << "returning value: " << best_val << " = " << dynamic_programming[best_index] << endl;
}

void PseudoTopoOrder::randomize()
{
	int start = 0;
	int end;
	for (const auto& X : m_parent.strongly_connected_components())
	{
		if (X.size() == 1)
			continue;
		end = start+X.size();
		random_shuffle(pto.begin()+start, pto.begin()+end);
		for (int i = start; i < end; ++i)
		{
			pto_inverse[pto[i]] = i;
		}
		start = end;
	}
	AnnounceModification(0);
}

void PseudoTopoOrder::transfer(int a, int b, int c, int d, int h)
{
	if (h == b-a)
	{
// 		cout << "PseudoTopoOrder::transfer nothing to do! abcdh = " << a << "," << b << "," << c << "," << d << ","  << h << endl;
		return;
	}
	
	while (h > b-a) // we must transfer from end to start
	{
		transpose(b,c);
		++b;
		++c;
	}
	
	
	while (h < b-a) // we must transfer from start to end
	{
		--b;
		--c;
		transpose(b,c);
	}
}

void PseudoTopoOrder::reverse_order(int a, int b)
{
	std::reverse(pto.begin()+a, pto.begin()+b);
	for (int i = a; i < b; ++i)
	{
		pto_inverse[pto[i]] = i;
	}
	AnnounceModification(a);
}

void PseudoTopoOrder::heuristic_sort(int a, int b, int numtimes)
{
	for (int r = b-1; r != a; --r)
	{
		node_t node = pto[r];
		int iu = get_exneighbor_in_range(a,r,node);
		if (iu != -1)
			transpose(iu,r);
	}
	
	for (int i = 0; i < numtimes; ++i)
	{
		int r = rand()%(b-a)+a;
		node_t node = pto[r];
		int iu = get_exneighbor_in_range(a,r,node);
		if (iu != -1)
			transpose(iu,r);
	}
}

int PseudoTopoOrder::get_exneighbor_in_range(int a, int b, node_t node)
{
	for (auto u : m_parent.exneighbors(node))
	{
		int iu = pto_inverse[u];
		if (a <= iu && iu < b)
			return iu;
	}
	return -1;
}

bool PseudoTopoOrder::eXtreme_edge_opener()
{
	auto oldval = Value();
	FillPath();

	for (const auto& x : m_parent.big_scc())
	{
		int a = x.start;
		int d = x.end;
		if (d - a < 5)
			continue;
		
		// true true true *false false
		auto f = std::partition_point(m_path.rbegin(), m_path.rend(), [this,a](int i) -> bool
		{
			return i < a;
		});

		int aa = a;
		while (*f < d && f != m_path.rend())
		{
			if (*f != aa)
				transpose(aa,*f);
			++f;
			++aa;
		}
		
		int c = d;
		int b = aa;
		
		if (c-b < 4)
			continue;
		
		auto order = range<int>(b-a+1);
		random_shuffle(order.begin(), order.end());
		shuffle(b,c);
		for (auto h : order)
		{
			transfer(a,b,c,d,h);
			int total = b-a + d-c;
			b = a+h;
			c = d-(total-h);
						
			heuristic_sort(b,c,15000);
			if (Value() > oldval)
			{
// 				cout << "(" << ChronometerPeek() << ", " << Value() << ")," << endl;
				FillPath();
				return true;
			}
		}

	}
	FillPath();
	return false;
}

void PseudoTopoOrder::open_edges_until_no_more_improvement_found(double maxnumseconds)
{
	clock_t start = clock();
	while (diffclock(clock(),start) < maxnumseconds)
	{
		eXtreme_edge_opener();
	}
}


void PseudoTopoOrder::FillPath()
{
	if (path_filled)
		return;
	m_path.clear();
	FillDP();
	auto m = Value();
	int a = best_index;
	
	while (true)
	{
		bool found = false;
// 		cout << "toret = " << toReturn << endl;
		node_t u = pto[a];
		m_path.push_back(a);
		for (auto v : m_parent.inneighbors(u))
		{
			node_t b = pto_inverse[v];
			if (b > a)
				continue;
			if (dynamic_programming[b] == dynamic_programming[a] - m_parent.edge_value(v,u))
			{
				a = b;
				m = dynamic_programming[b];
				found = true;
				break;
			}
		}
		if (!found)
		{
			path_filled = true;
			return;
		}
	}
}
