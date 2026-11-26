# Longest Simple Path

This attempts to find a long simple path (without repeated vertices) in a weighted digraph. It uses heuristics to find good paths.

# Quick preview:

```c++
    DiGraph D(10000); // creates a digraph on 10000 vertices named 0, 1, ... , 9999
    D.add_edge(1,27,3); // adds an edge to the digraph from node 1 to node 27 with weight 3
    Path P = D.FindLongestSimplePath(); // Uses 5 seconds to find a long simple path
```
prints out the path.

Vertices can be named too:
    
```c++
    DiGraph G({"A","B","C","D","E","F","G","H"});
    
    G.add_edge("A","B",5);
    G.add_edge("B","C",2);
    G.add_edge("C","A",3);
    G.add_edge("C","D",4);
    G.add_edge("D","E",1);
    G.add_edge("B","F",3);
    G.add_edge("F","G",2);
    
    cout << G << endl;
    Path PG = G.FindLongestSimplePath();
    
    cout << "The best path I found has value " << 
PG.Value() << endl;
```

Of course, this program is not optimized at all for small graphs, but for large (~10,000 vertices) graphs. It will consume the time it has trying to improve something that can't be improved...

# Note 1
Graph D gets modified when calling FindLongestSimplePath, so beware. The problem is that we erase nodes in "small" connected components. If possible, make sure D is (weakly) connected. Else, the program will only focus on the largest (weakly) connected component. For example, is you have two components: A->B->C and D->E, but weight of D->E is 1000 and weight of A->B and B->C are both 1, it will ignore the D->E, since it's a small component. For most graphs this doesn't matter, since the longest path will likely be in the largest connected component.

# Improving your results

It's very likely you can improve your results. Many choices were made in order to improve stability of the algorithm (i.e. repeat some searches a few times, and so on) that can be removed and in most cases work better, but are easier to "trick".
See class SearchOptions for a bunch of options.

# The method

The method it uses will be described in an upcoming research paper (in process), but it does a lot of things. You can also invite me to give a talk on this.

# Acknowledgements

Research supported in part by PAPIIT IA106316, UNAM.


