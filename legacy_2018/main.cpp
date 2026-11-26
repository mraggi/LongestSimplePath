#include "digraph.hpp"
#include "digraph.hpp"
#include "path.hpp"
#include "trainer.hpp"
#include <cmath>

int main() 
{
	std::ios_base::sync_with_stdio(false);
	srand(time(nullptr));
// 	randomize();
	TimeFromStart(); // in order to create the timer. Doesn't matter much.
	// To create a graph with names, just pass a vector<string> containing names to the constructor.
    DiGraph G({"A","B","C","D","E","F","G","H"});
    
    G.add_edge("A","B",7);
    G.add_edge("B","C",2);
    G.add_edge("C","A",3);
    G.add_edge("C","D",4);
    G.add_edge("D","E",4);
    G.add_edge("B","F",3);
    G.add_edge("F","G",2);

	std::cout << "Created graph!" << std::endl;
    std::cout << G << std::endl;
    Path PG = G.FindLongestSimplePath();
    
    std::cout << "The best path I found for graph G with the default options has value " << PG.value() << std::endl;

    
    int n = 2000;
    double p = log(double(n))/(1.0*n);
    
    std::cout << "\n\nCreating Erdos Renyi graph..." << std::flush;
    // Erdos Renyi digraph with n vertices, where each edge has a probability of p of appearing.
    DiGraph D = DiGraph::CreateRandomDiGraph(n,p);
    std::cout << " done!" << "\nSearching Path..." << std::endl;
    
    
    Path P = D.FindLongestSimplePath();
    
    std::cout << "The best path I found for graph D has value " << P.value() << '\n';
	
    return 0;
}
