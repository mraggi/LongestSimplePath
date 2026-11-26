#pragma once
#include <vector>
#include <unordered_map>
#include <map>
#include <iostream>
#include <sstream>
#include <chrono>
#include <ctime>
#include <random>
#include <fstream>
#include <cassert>

template <class T> 
T max(const std::vector<T>& v)
{
	assert (v.size() != 0);
	T m = v[0];
	size_t sz = v.size();
	for (size_t i = 0; i < sz; ++i)
	{
		if (v[i] > m)
			m = v[i];
	}
	return m;
}

template <class T> std::vector<T> range (T n) 
{
	std::vector<T> toReturn(n, 0);
	for (T i = 1; i < n ; ++i)
		toReturn[i] = i;
	return toReturn;
}

template <class T> 
	std::ostream& operator<<(std::ostream& os, const std::vector<T>& rhs)
	{
		os << "[ ";
		for (const auto& x : rhs)
			os << x << " ";
		os << "]";
		return os;
	}
	
	/*
template <class T>
class SquareSparseMatrix
{
public:
    SquareSparseMatrix(int i) : n(i)
	{ 
		m_data = std::vector<std::vector<T>>(n,std::vector<T>(n,0));
	}
    
//     inline int size() const { return n; }
    
    inline     int to2d(int x, int y) const { return x*n+y; }
    
    inline
    T& operator()(int x, int y)
    {
// 		return m_data[to2d(x,y)];
        return m_data[x][y];
    }
    
    inline
    T operator()(int x, int y) const
    {
// 		auto f = m_data.find(to2d(x,y));
//         return f->second;
// 		return m_data[to2d(x,y)];
		return m_data[x][y];
    }
    
//     inline
//     T operator()(size_t x, size_t y) const
//     {
//         return m_data[to2d(x,y)];
//     }

private:
//     std::vector<std::vector<T>> m_data;
	std::vector<std::vector<T>> m_data;
    int n;
};

template <class T> 
std::ostream& operator<<(std::ostream& os, const SquareSparseMatrix<T>& M)
{
	for (size_t x = 0; x < M.size(); ++x)
	{
		for (size_t y = 0; y < M.size(); ++y)
		{
			os << M(x,y) << " ";
		}
		os << std::endl;
	}
	return os;
}*/

// inline double diffclock(clock_t a, clock_t b)
// {
// 	const double c = 1.0/CLOCKS_PER_SEC;
// 	return double(a-b)*c;
// }

typedef std::chrono::time_point<std::chrono::high_resolution_clock> clockt;

inline double diffclockt(clockt a, clockt b)
{
	
	const double t = 0.000001;
	return std::chrono::duration_cast<std::chrono::microseconds>(a-b).count()*t;
}

class Chronometer
{
public:
	Chronometer() : m_timer(std::chrono::high_resolution_clock::now()) {}
	
	double Reset()
	{
		auto tlast = m_timer;
		m_timer = std::chrono::high_resolution_clock::now();
	
		return diffclockt(m_timer, tlast);
	}
	
	double Peek() const
	{
		auto tnow = std::chrono::high_resolution_clock::now();
	
		return diffclockt(tnow, m_timer);
	}
	
	std::chrono::time_point<std::chrono::high_resolution_clock> m_timer;
};

inline double TimeFromStart()
{
	static Chronometer C;
	return C.Peek();
}

/*
inline double Chronometer()
{
	auto tlast = RClock::Instance().running_timer;
	RClock::Instance().running_timer = std::chrono::high_resolution_clock::now();
	return diffclockt(RClock::Instance().running_timer,tlast);
}

inline double ChronometerPeek()
{
	return diffclockt(std::chrono::high_resolution_clock::now(),RClock::Instance().running_timer);
}*/

inline std::default_random_engine & random_engine()
{
	static std::default_random_engine e{};
	return e;
}

inline void randomize()
{
	static std::random_device rd{};
	random_engine().seed(rd());
}

inline bool probability_of_true(double p)
{
	static std::bernoulli_distribution d(p);
	return d(random_engine());
}

inline double random_real( double from, double upto )
{
	static std::uniform_real_distribution<> d{};
	using parm_t = decltype(d)::param_type;
	return d( random_engine(), parm_t{from, upto} );
}

inline int random_integer(int from, int to)
{
	static std::uniform_int_distribution<> d{};
	using parm_t = decltype(d)::param_type;
	return d( random_engine(), parm_t{from, to-1} );
}

inline int random_give_priority_to_primeros(int a, int b)
{
	int n = b-a;
	random_engine();
	int t = random_integer(0,(n*(n+1))/2);
	int u = n-1;
	int i = n-1;
	int toreturn = a;
	while (u < t && i > 0)
	{
		u += i;
		--i;
		++toreturn;
	}
	return toreturn;
}

std::vector<std::string> openFile(const std::string& filename);

template<class Predicate>
std::vector<std::string> openFile(const std::string& filename, Predicate pred)
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
		if (pred(line))
			toreturn.push_back(line);
    }
    
//     clock_t read = clock();
// 	std::cout << "Time to read file: " << diffclock(read,start) << std::endl;
	return toreturn;
}

std::vector<std::string> readfromstdin();
