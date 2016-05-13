#include "digraph.hpp"
#include "path.hpp"
#include "moviegraph.hpp"
#include "trainer.hpp"
#include <cmath>

using namespace std;

int main() 
{
	srand(time(nullptr));
	randomize();
	TimeFromStart(); // in order to create the timer. Doesn't matter much.
	int N = 8;
	
//     auto V = openFile("all_titles.txt");
//     
//     auto M = new MovieGraph(V);
//     M->process();
//     
// 	cout << "V(M) = " << M->num_vertices() << " and E(M) = " << M->num_edges() << endl;
// 	int VM = M->num_vertices();
// 	double avgedgeweight = 0;
// 	for (int i = 0; i < VM; ++i)
// 	{
// 		for (int j = 0; j < VM; ++j)
// 		{
// 			avgedgeweight+=M->edge_value(i,j);
// 		}
// 	}
// 	cout << "avg edge weight: " << avgedgeweight/M->num_edges() << endl;
	
// 	return 0;
// 	cout << "average edge weight: " << 
	
	const int n = 3000;
	const double p = (log(double(n)))/(n);
// 	const double p = double(M.num_edges())/(n*(n-1));
	vector<DiGraph*> graphs;
// 	graphs.push_back(M);
	
	for (int i = 0; i < 30; ++i)
	{
		cout << endl << "Creating graph " << i << endl;
		auto D = new DiGraph(DiGraph::CreateRandomWeightedDiGraph(n,p,3,20));
		D->process();
		cout << "V(D) = " << D->num_vertices() << " and E(D) = " << D->num_edges() << endl;
		graphs.push_back(D);
		D->FindLongestSimplePath(1.0);
	}
	
	Trainer T(graphs);
	
	auto a = T.Train(15);
	cout << "after training, I get these parameters: " << a << endl;
	
	double totalImprovementA = 0;
	double totalImprovementB = 0;
	double totalImprovementC = 0;
	
    int numcrossvalidation = 1000;
    
	for (int i = 0; i < numcrossvalidation; ++i)
	{
		cout << endl << "Cross-validating: " << i << endl;
		DiGraph D(DiGraph::CreateRandomWeightedDiGraph(n,p,3,20));
		
		auto beforeP = D.FindLongestSimplePath(1.0).Value();
		cout << "cross-validation BEFORE " << beforeP << endl;
		D.set_parameters(a);
		auto afterP = D.FindLongestSimplePath(1.0).Value();
		cout << "cross-validation AFTER " << afterP << endl;
		D.set_parameters({1,0,0,0,1,0,0,0});
		auto controlP = D.FindLongestSimplePath(1.0).Value();
		cout << "Control 00000000: " << controlP << endl;
		D.set_parameters({1,1,1,1,1,1,1,1});
		auto control2P = D.FindLongestSimplePath(1.0).Value();
		cout << "Control 11111111: " << control2P << endl;
		cout << "\tImprovement over 1 4 16 64: " << afterP-beforeP << endl;
		cout << "\tImprovement over 1 0 0 0: " << afterP-controlP << endl;
		cout << "\tImprovement over 1 1 1 1: " << afterP-control2P << endl;
		totalImprovementA += afterP-beforeP;
		totalImprovementB += afterP-controlP;
		totalImprovementC += afterP-control2P;
	}
	
	cout << "DONE! Total improvement over  1 4 16 64: " << totalImprovementA/numcrossvalidation << endl;
	cout << "DONE! Total improvement over  1 0 0 0: " << totalImprovementB/numcrossvalidation << endl;
	cout << "DONE! Total improvement over  1 1 1 1: " << totalImprovementC/numcrossvalidation << endl;
	
	for (auto p : graphs)
		delete p;
	
    return 0;
}
