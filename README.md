# Longest Simple Path

This attempts to find a long simple path (without repeated vertices) in a weighted digraph. It uses heuristics to find good paths.

# Quick preview:

```c++
    DiGraph D(10000); // creates a digraph on 10000 vertices named 0, 1, ... , 9999
    D.add_edge(1,27,3); // adds an edge to the digraph from node 1 to node 27 with weight 3
    Path P = D.FindLongestSimplePath(5.0); // Uses 5 seconds to find a long simple path
    cout << P << endl;
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
    Path PG = G.FindLongestSimplePath(8.0);
    
    cout << "The best path I found in one second has value " << PG.Value() << " and is " << PG << endl;
```

# Note
Graph D gets modified when calling FindLongestSimplePath, so beware. The problem is that we erase nodes in "small" connected components. If possible, make sure D is (weakly) connected. Else, the program will only focus on the largest (weakly) connected component. For example, is you have two components: A->B->C and D->E, but weight of D->E is 1000 and weight of A->B and B->C are both 1, it will ignore the D->E, since it's a small component. For most graphs this doesn't matter, since the longest path will likely be in the largest connected component.

# The method

The method it uses will be described in an upcoming research paper (in process), but it does a lot of things. You can also invite me to give a talk on this.