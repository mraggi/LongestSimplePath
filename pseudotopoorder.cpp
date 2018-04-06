#include "pseudotopoorder.hpp"


void PseudoTopoOrder::shuffle(int a, int b)
{
	random_shuffle(pto.begin() + a, pto.begin() + b);

	for (int i = a; i < b; ++i)
	{
		pto_inverse[pto[i]] = i;
	}

	AnnounceModification(a);
}

void PseudoTopoOrder::RecalcTopoInverse()
{
	for (int i = 0; i < pto.size(); ++i)
	{
		pto_inverse[pto[i]] = i;
	}

	AnnounceModification(0);
}

void PseudoTopoOrder::apply(const Path& P)
{
	global_best = P.value();
	auto it = P.get_path().begin();
	int fu = -1;
	int n = pto.size();

	for (int i = 0; i < n; ++i)
	{
		int x = pto[i];

		if (P.is_node_in_path(x))
		{
			if (fu == -1)
			{
				fu = i;
			}

			pto[i] = *it;
			pto_inverse[*it] = i;
			++it;
		}
	}
	AnnounceModification(fu);
// 	FillDP();
// 	FillPath();
	
}

void PseudoTopoOrder::apply(const Path& P, int u, int v)
{
// 	cout << "Before applying Q, we get: " << Value() << endl;

	std::vector<int> indexesofPbetweenuandv;
	std::vector<int> nodesofPbetweenuandv;

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
// 	std::cout << "Filling path, with value: " << Value() << std::endl;
	Path P(size());

	for (auto i : m_path)
	{
		P.emplace_front(pto[i], i.Weight());
	}

	return P;
}

void PseudoTopoOrder::FillDP()
{
	int n = pto.size();
	int best_val = 0;

	if (best_index < first_unknown)
	{
		best_val = dynamic_programming[best_index];
	}

// 	std::cout << "\nAt the beginning of filldp, first unknown: " << first_unknown << "/" << n << " and value of bestindex = " << dynamic_programming[best_index] << std::endl;
	for ( ; first_unknown < n; ++first_unknown)
	{
		int u = pto[first_unknown];
		dynamic_programming[first_unknown] = 0;

		auto& neigh = m_parent.inneighbors(u);

		for (auto v : neigh)
		{
			auto j = pto_inverse[v];

			if (first_unknown <= j)   // should be ignored
			{
				continue;
			}

			auto candidate = dynamic_programming[j] + v.Weight();

			if (candidate > dynamic_programming[first_unknown])
			{
				dynamic_programming[first_unknown] = candidate;
				best_parent[first_unknown] = j;

				if (candidate > best_val)
				{
// 					cout << "Fill: new best value at " << i << " with value " << best_val << endl;
					best_val = candidate;
					best_index = first_unknown;
				}
			}
		}
	}
	if (best_val > global_best)
	{
		std::cout << "PTO improved path to " << best_val << " at " << TimeFromStart() << std::endl;
		global_best = best_val;
	}

// 	std::cout << "\nreturning value: " << best_val << " = " << dynamic_programming[best_index] << std::endl;
}

void PseudoTopoOrder::randomize()
{
	int start = 0;
	int end;

	for (const auto& X : m_parent.strongly_connected_components())
	{
		if (X.size() == 1)
		{
			continue;
		}

		end = start + X.size();
		random_shuffle(pto.begin() + start, pto.begin() + end);

		for (int i = start; i < end; ++i)
		{
			pto_inverse[pto[i]] = i;
		}

		start = end;
	}

	AnnounceModification(0);
}

void PseudoTopoOrder::transfer(int a, int b, int c, int  /*d*/, int h)
{
	if (h == b - a)
	{
// 		cout << "PseudoTopoOrder::transfer nothing to do! abcdh = " << a << "," << b << "," << c << "," << d << ","  << h << endl;
		return;
	}

	while (h > b - a) // we must transfer from end to start
	{
		transpose(b, c);
		++b;
		++c;
	}


	while (h < b - a) // we must transfer from start to end
	{
		--b;
		--c;
		transpose(b, c);
	}
}

void PseudoTopoOrder::reverse_order(int a, int b)
{
	std::reverse(pto.begin() + a, pto.begin() + b);

	for (int i = a; i < b; ++i)
	{
		pto_inverse[pto[i]] = i;
	}

	AnnounceModification(a);
}

void PseudoTopoOrder::heuristic_sort(int a, int b, int numtimes)
{
	for (int r = b - 1; r != a; --r)
	{
		node_t node = pto[r];
		int iu = get_outneighbor_in_range(a, r, node);

		if (iu != -1)
		{
			transpose(iu, r);
		}
	}

	for (int i = 0; i < numtimes; ++i)
	{
		int r = rand() % (b - a) + a;
		node_t node = pto[r];
		int iu = get_outneighbor_in_range(a, r, node);

		if (iu != -1)
		{
			transpose(iu, r);
		}
	}
}

int PseudoTopoOrder::get_outneighbor_in_range(int a, int b, node_t node)
{
	for (auto u : m_parent.outneighbors(node))
	{
		int iu = pto_inverse[u];

		if (a <= iu && iu < b)
		{
			return iu;
		}
	}

	return -1;
}

bool PseudoTopoOrder::is_path_in_order(const Path& P)
{
	auto it = P.get_path().begin();
	auto it2 = it;
	++it2;
	while ( it2 != P.get_path().end())
	{
		if (pto_inverse[*it] >= pto_inverse[*it2])
			return false;
		++it;
		++it2;
		
	}
	return true;
}


bool PseudoTopoOrder::eXtreme_edge_opener()
{
	auto oldval = Value();
	FillPath();
	int max_pointless = m_parent.Options.pto_scc_size_max_pointless;
	for (const auto& x : m_parent.big_scc())
	{
		int a = x.start;
		int d = x.end;

		if (d - a <= max_pointless)
		{
			continue;
		}

		// true true true *false false
		auto f = std::partition_point(m_path.rbegin(), m_path.rend(), [this, a](node_t i) -> bool
		{
			return i < a;
		});

		int aa = a;

		while (*f < d && f != m_path.rend())
		{
			if (*f != aa)
			{
				transpose(aa, *f);
			}

			++f;
			++aa;
		}

		int c = d;
		int b = aa;

		if (c - b < max_pointless)
		{
			continue;
		}

		auto order = range<int>(b - a + 1);
		random_shuffle(order.begin(), order.end());
		shuffle(b, c);
		int ntimeshsort = m_parent.Options.pto_num_heuristic_sort;
		for (auto h : order)
		{
			transfer(a, b, c, d, h);
			int total = b - a + d - c;
			b = a + h;
			c = d - (total - h);

			heuristic_sort(b, c, ntimeshsort);

			if (Value() > oldval)
			{
				timer.Reset();
				FillPath();
				return true;
			}
		}

	}

	FillPath();
	return false;
}

void PseudoTopoOrder::open_edges_until_no_more_improvement_found()
{
	double maxnumseconds = m_parent.Options.pto_time_without_improvement;

	while (timer.Peek() < maxnumseconds)
	{
		eXtreme_edge_opener();
	}
}


void PseudoTopoOrder::FillPath()
{
	if (path_filled)
	{
		return;
	}

	m_path.clear();
	FillDP();
// 	auto m = Value();
// 	std::cout << "\nAt beginning of fillpath, value: " << Value() << std::endl;
	weight_t currweight = 0;
	int a = best_index;
	path_filled = true;
	int numadded = 0;
	while (true)
	{
		bool found = false;
// 		cout << "toret = " << toReturn << endl;
		node_t u = pto[a];
		m_path.emplace_back(a, currweight);
		++numadded;
// 		std::cout << "Adding " << a << " to path" << std::endl;
		for (auto v : m_parent.inneighbors(u))
		{
			node_t b = pto_inverse[v];

// 			if (b > a)
// 			{
// // 				continue;
// 			}

			if (b < a && dynamic_programming[b] + v.Weight() == dynamic_programming[a])
			{
				a = b;
// 				m = dynamic_programming[b];
				currweight = v.Weight();
				found = true;
				break;
			}
		}

		if (!found)
		{
// 			std::cout << "\nAt end of fillpath, value: " << numadded << std::endl;
			if (dynamic_programming[a] != 0)
			{
				std::cout << "\n\n\n AAAAAAAAAAAAAAAAAA: I still have " << dynamic_programming[a] << " to go at vertex at index " << a << " but I can't find where to go from here! I'm supposed to go to " << best_parent[a] << std::endl;
				node_t u = pto[a];
				for (auto v : m_parent.inneighbors(u))
				{
					auto b = pto_inverse[v];
					if (b < a)
					{
						std::cout << "\t" << b << " with value " << dynamic_programming[b] << std::endl;
					}
				}
				throw;
			}
			return;
		}
	}
}
