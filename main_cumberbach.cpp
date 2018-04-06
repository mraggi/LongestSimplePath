#include "digraph.hpp"
#include "digraph.hpp"
#include "path.hpp"
#include "trainer.hpp"
#include <cmath>

using namespace std;

int main() 
{
	int n = 74;
	DiGraph G(n);
// 	vector<vector<int>> A(n,vector<int>(n,0));
	for (int x = 0; x < n; ++x)
	{
		for (int y = 0; y < n; ++y)
		{
// 			cin >> A[x][y];
			int a;
			cin >> a;
			if (a != 0)
				G.add_edge(x,y,a);
		}
	}
	
// 	cout << "Terminé de leer, con " <<  << endl;
	
	Path P = G.FindLongestSimplePath();
	
	cout << P.value() << endl;
	
	
    return 0;
}
