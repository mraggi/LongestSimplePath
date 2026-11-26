#include "helpers.hpp"
#include <fstream>

std::vector<std::string> openFile(const std::string& filename)
{
	std::vector<std::string> toreturn;
	std::string line;
    std::ifstream myfile(filename);

    if(!myfile) //Always test the file open.
    {
        std::cout<<"Could not open file " << filename << std::endl;
		throw "File not found";
    }
    
//     clock_t start = clock();  
	
    while (std::getline(myfile, line))
    {
        toreturn.push_back(line);
    }
    
//     clock_t read = clock();
// 	std::cout << "Time to read file: " << diffclock(read,start) << std::endl;
	return toreturn;
}

std::vector<std::string> readfromstdin()
{
// 	clock_t start = clock();
	std::string bla;
	std::vector<std::string> toreturn;
//  	toreturn.reserve(10000);
	while (std::cin)
	{
		getline(std::cin,bla);
		toreturn.push_back(bla);
	}
	return toreturn;
}
