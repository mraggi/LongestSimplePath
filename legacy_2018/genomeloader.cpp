#include "digraph.hpp"
#include "path.hpp"
#include "genomegraph.hpp"
#include "trainer.hpp"
#include <cmath>

using namespace std;

int main() 
{
	srand(time(nullptr));
	randomize();
	TimeFromStart(); // in order to create the timer. Doesn't matter much.
	
    auto f = openFile("fragmentos.fasta", [](const string& s){
		return s[0] != '>';
	});
    sort(f.begin(), f.end());
// 	for (auto t : f) cout << t << endl;
// 	return 0;
	cout << "creating graph..." << endl;
    GenomeGraph M(f);
	cout << "Done loading graph! Num vertices: " << M.num_vertices() << " and edges = " << M.num_edges() << endl;
	
// 	return 0;
    
    Path P = M.FindLongestSimplePath(1.0);
	
    cout << "Value = " << P.value() << endl;
// 	M.PrintPath(P);
	cout << "Time taken: " << TimeFromStart() << endl;
//     cout << "Path = " << P << endl;
    return 0;
}
