#include "helpers.hpp"
#include <fstream>

vector<string> openFile(const string& filename)
{
	vector<string> toreturn;
	string line;
    std::ifstream myfile(filename);

    if(!myfile) //Always test the file open.
    {
        cout<<"Could not open file " << filename << endl;
		throw "File not found";
    }
    
    clock_t start = clock();
	
    while (std::getline(myfile, line))
    {
        toreturn.push_back(line);
    }
    
    clock_t read = clock();
// 	cout << "Time to read file: " << diffclock(read,start) << endl;
	return toreturn;
}

vector<string> readfromstdin()
{
	clock_t start = clock();
	string bla = "";
	vector<string> toreturn;
//  	toreturn.reserve(10000);
	while (cin)
	{
		getline(cin,bla);
		toreturn.push_back(bla);
	}
	return toreturn;
}