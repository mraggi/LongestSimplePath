# Longest Simple Path

This attempts to find a long simple path (without repeated vertices) in a weighted digraph. It uses heuristics to find good paths.

# Quick preview:

```c++
    DiGraph D(10000); // creates a digraph on 10000 vertices named 0, 1, ... , 9999
    D.add_edge(1,27,0.8); // adds an edge to the digraph from node 1 to node 27 with weight 0.8
    Path P = D.FindLongestSimplePath(5.0); // Uses 5 seconds to find a long simple path
	auto dq = P.get_path();
```

dq will be a `deque<int>` that contains the nodes. `P.Value()` returns the total weight.

# Note
Graph D gets modified when calling FindLongestSimplePath, so beware. The returned path contains the old nodes, but don't try to call use D for other purposes.

# The method

The method it uses will be described in an upcoming research paper (in process), but it does a lot of things. The file "talk.pdf" contains a talk in spanish that somewhat explains it.

