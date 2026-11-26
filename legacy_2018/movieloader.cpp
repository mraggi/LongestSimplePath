#include "digraph.hpp"
#include "path.hpp"
#include "moviegraph.hpp"
#include "trainer.hpp"
#include <cmath>

int main() 
{
	srand(time(nullptr));
	randomize();
	TimeFromStart(); // in order to create the timer. Doesn't matter much.
	
    auto f = openFile("all_titles.txt");
    sort(f.begin(), f.end());
    MovieGraph M(f);
    
    Path P = M.FindLongestSimplePath(1.0);
    cout << "Value = " << P.value()/2 << endl;
	M.PrintPath(P);
	cout << "Time taken: " << TimeFromStart() << endl;
//     cout << "Path = " << P << endl;
    return 0;
}
