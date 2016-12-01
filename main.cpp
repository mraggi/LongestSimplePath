#include "digraph.hpp"
#include "path.hpp"
#include "moviegraph.hpp"
#include "trainer.hpp"
#include <cmath>

using namespace std;

int main() 
{
// 	srand(time(nullptr));
// 	randomize();
	TimeFromStart(); // in order to create the timer. Doesn't matter much.
	
    DiGraph G({"A","B","C","D","E","F","G","H"});
    
    G.add_edge("A","B",7);
    G.add_edge("B","C",2);
    G.add_edge("C","A",3);
    G.add_edge("C","D",4);
    G.add_edge("D","E",4);
    G.add_edge("B","F",3);
    G.add_edge("F","G",2);
    
    cout << G << endl;
    Path PG = G.FindLongestSimplePath(0.5);
    
    cout << "The best path I found in half a second has value " << PG.value() << endl;

    
    int n = 1001;
    double p = log(double(n))/(1.0*n);
    
    cout << "Creating graph..." << flush;
    // Erdos Renyi digraph with n vertices, where each edge has a probability of p of appearing.
    DiGraph D = DiGraph::CreateRandomDiGraph(n,p);
    cout << " done!" << endl;
    
    cout << "Processing graph...";
    D.process();
    cout << " done!" << endl;
    
    cout << "Searching path...";    
    Path P = D.FindLongestSimplePath(1.0,0.05,4);
    cout << " done!" << endl;
    
    cout << "The best path I found has value " << P.value() << endl;
	
    return 0;
}
