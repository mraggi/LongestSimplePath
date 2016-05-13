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
	
    auto f = openFile("all_titles.txt");
    
    MovieGraph M(f);
    
    Path P = M.FindLongestSimplePath(1.0);
    cout << "Value = " << P.Value() << endl;
    cout << "Path = " << P << endl;
    return 0;
}
