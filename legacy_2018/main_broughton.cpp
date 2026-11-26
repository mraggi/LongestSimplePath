#include "digraph.hpp"
#include "digraph.hpp"
#include "path.hpp"
#include "trainer.hpp"
#include <cmath>

using namespace std;

int main() 
{
	int n = 135;
	DiGraph G(n);
// 	vector<vector<int>> A(n,vector<int>(n,0));
	G.add_edge(1,12, 244);
G.add_edge(2,6, 214);
G.add_edge(2,3, 10);
G.add_edge(2,10, 33);
G.add_edge(4,11, 41);
G.add_edge(5,14, 78);
G.add_edge(7,8, 14);
G.add_edge(8,9, 32);
G.add_edge(8,15, 74);
G.add_edge(9,10, 61);
G.add_edge(9,20, 34);
G.add_edge(10,12, 38);
G.add_edge(11,13, 4);
G.add_edge(11,14, 29);
G.add_edge(12,16, 9);
G.add_edge(14,15, 17);
G.add_edge(15,19, 42);
G.add_edge(16,17, 8);
G.add_edge(16,24, 72);
G.add_edge(18,21, 3);
G.add_edge(20,21, 20);
G.add_edge(20,25, 34);
G.add_edge(21,23, 20);
G.add_edge(22,26, 30);
G.add_edge(25,27, 62);
G.add_edge(25,30, 23);
G.add_edge(26,31, 24);
G.add_edge(26,36, 85);
G.add_edge(28,31, 8);
G.add_edge(29,30, 53);
G.add_edge(30,35, 31);
G.add_edge(31,36, 67);
G.add_edge(32,33, 15);
G.add_edge(33,34, 31);
G.add_edge(33,35, 31);
G.add_edge(35,39, 88);
G.add_edge(36,38, 48);
G.add_edge(37,40, 65);
G.add_edge(38,39, 40);
G.add_edge(38,41, 14);
G.add_edge(39,40, 2);
G.add_edge(40,43, 23);
G.add_edge(41,42, 64);
G.add_edge(41,45, 21);
G.add_edge(43,46, 64);
G.add_edge(43,53, 48);
G.add_edge(44,49, 33);
G.add_edge(47,51, 96);
G.add_edge(48,51, 70);
G.add_edge(49,52, 4);
G.add_edge(49,56, 34);
G.add_edge(50,55, 12);
G.add_edge(51,52, 30);
G.add_edge(52,62, 51);
G.add_edge(53,54, 3);
G.add_edge(53,63, 99);
G.add_edge(54,55, 28);
G.add_edge(54,57, 10);
G.add_edge(55,59, 14);
G.add_edge(58,67, 73);
G.add_edge(60,67, 99);
G.add_edge(61,65, 54);
G.add_edge(62,64, 39);
G.add_edge(62,66, 18);
G.add_edge(63,68, 44);
G.add_edge(63,91, 105);
G.add_edge(65,66, 58);
G.add_edge(65,71, 41);
G.add_edge(66,72, 41);
G.add_edge(67,70, 5);
G.add_edge(69,70, 21);
G.add_edge(70,73, 14);
G.add_edge(72,76, 59);
G.add_edge(72,80, 40);
G.add_edge(73,74, 17);
G.add_edge(73,81, 37);
G.add_edge(75,79, 24);
G.add_edge(77,78, 5);
G.add_edge(78,79, 34);
G.add_edge(78,82, 17);
G.add_edge(79,81, 27);
G.add_edge(80,88, 37);
G.add_edge(80,96, 48);
G.add_edge(81,93, 16);
G.add_edge(83,99, 68);
G.add_edge(84,85, 42);
G.add_edge(85,86, 56);
G.add_edge(85,97, 38);
G.add_edge(87,98, 83);
G.add_edge(89,90, 18);
G.add_edge(90,91, 7);
G.add_edge(90,100, 38);
G.add_edge(91,92, 33);
G.add_edge(93,94, 33);
G.add_edge(93,103, 39);
G.add_edge(95,104, 16);
G.add_edge(96,98, 2);
G.add_edge(96,99, 143);
G.add_edge(97,100, 135);
G.add_edge(97,101, 28);
G.add_edge(98,110, 45);
G.add_edge(99,101, 51);
G.add_edge(100,102, 25);
G.add_edge(101,113, 51);
G.add_edge(102,103, 63);
G.add_edge(102,111, 45);
G.add_edge(103,105, 6);
G.add_edge(104,105, 89);
G.add_edge(104,107, 10);
G.add_edge(105,121, 56);
G.add_edge(106,111, 30);
G.add_edge(108,110, 37);
G.add_edge(109,123, 80);
G.add_edge(110,118, 44);
G.add_edge(111,119, 107);
G.add_edge(112,113, 107);
G.add_edge(113,114, 1);
G.add_edge(114,115, 108);
G.add_edge(114,116, 21);
G.add_edge(116,122, 134);
G.add_edge(116,125, 39);
G.add_edge(117,118, 36);
G.add_edge(118,124, 22);
G.add_edge(120,121, 12);
G.add_edge(121,126, 16);
G.add_edge(123,124, 46);
G.add_edge(123,129, 31);
G.add_edge(124,128, 27);
G.add_edge(127,128, 36);
G.add_edge(128,130, 24);
G.add_edge(130,131, 37);
G.add_edge(130,133, 36);
G.add_edge(132,133, 25);
G.add_edge(133,134, 2);
	
// 	cout << "Terminé de leer, con " <<  << endl;
	
	Path P = G.FindLongestSimplePath();
	
	cout << P.value() << endl;
	
	
    return 0;
}
