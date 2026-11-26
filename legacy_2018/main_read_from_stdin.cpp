#include "digraph.hpp"
#include "digraph.hpp"
#include "path.hpp"
#include "trainer.hpp"
#include <cmath>

DiGraph ReadDigraphFromStdin()
{
	using std::cin;
	
	int n,m;
	cin >> n >> m;
	DiGraph D(n);
	for (int i = 0; i < m; ++i)
	{
		int a, b;
		cin >> a >> b;
		D.add_edge(a,b);
	}
	return D;
}

int main() 
{
	using namespace std;
	std::ios_base::sync_with_stdio(false);
	
	auto D = ReadDigraphFromStdin();
	
	D.FindLongestSimplePath();
	
    return 0;
}
