Here’s a very brief, simple, and intuitive explanation of the solution to algo challange 10 
that refers back to original paper for details and proofs. 

We need 2 paths of different lengths. One path is the shortest path from s to t. 
We find this shortest path by doing BFS from s to t. To get the longer path, 
we first need to find all the shortest paths from s to t: 

Lets find all the shortest paths from s to t in the graph G. 
The edges and vertices of all the shortest paths from s to t 
will make up subgraph K. Subgraph K can be constructed in polynomial time 
(look at original paper or comments in polynomial_solution.py for create_shortest_paths_dag method. 

It does this using 2 BFS’s and for looping through every vertex.). 
Subgraph K is also a directed acyclic graph from S to T. 

In the code we call subgraph K the shortest_paths_dag. 
2 vertices can be in K but an edge between those 2 vertices might be missing
even though that edge exists in the original graph. That edge is not part of any shortest paths between 
s and t, so it isn't in shortest_paths_dag. However we refer to these edges as "LOST" edges because they didn't get added to K.
LOST edges can be used to construct longer paths. 


There are 2 ways to construct a longer path. The lost edges method and the outer vertex method. 
Both algorithms need to be used to say without doubt that there is no longer path in a graph. 
The lost edges method is described in the original paper. It tries to create a longer path 
using only the vertices in K. If we can find an edge E between 2 vertices in K, 
and edge E is not in the edge set of K, then E is a “lost edge” in K. E is not part of 
any shortest path, that is why it was not added to K, so E cannot be a 
forward edge or a tree edge. A simple path that traverses the lost edge to go 
from S to T is a longer path (proof, pseudocode, and detailed explanation in the original paper. In the code,
in polynomial_solution.py, lost edges method is 
implemented as the method called create_longer_path_using_lost_edges ).  

Since the lost edge is a back edge, we are going backward in K, and this might 
fail if no simple path can be made from s to t as a result of going backward. 
The method in the code:
    merge_two_overlapping_paths_in_dag(s, t, X, Y, shortest_paths_dag)
checks if using the backedge is ok by merging the paths s->V and F->t (the lost edge is V->F) 
using only vertices in the shortest paths dag. If it returns 2 simple paths, then we have found a longer path!
This method does a bfs to create s->V and then a dfs to create F->t and if that fails, tries to do a bfs from F->t and then a 
dfs from s->V and if that fails, returns false for being able to find those 2 simple paths in the DAG.

This method is also used in outer vertex method
(the second way to find a longer path).



    
The second way to find a longer path is using the outer vertex method. 
To find a longer path, we need to find a vertex F that isn't in a shortest path  
and attempt to construct a longer path by going from S to F, and then F to T. 
In other words, to find a longer path, we need to identify a vertex F called 
the coordination point in G-K that we can use so that we can go from S to F, 
then F to T (we denote this longer path [S,F,T]). This will be a longer path 
because it cannot be a shortest path; otherwise, the edges 
and vertices of [S,F,T] are already in K.


We can test all the vertices in G-K to see if any of them are coordination
points for valid longer simple paths. Creating the longer path with a 
coordination point is explained thoroughly in detail with proofs in the original paper, 
and can be done in polynomial time. 

Here is an intuitive explanation for the check:
Since we want a path from S to F to T, we will need to use BFS somehow to 
identify the paths from S to F, and a path from F to T, and then we join 
them to make the path [S,F,T]. Since these paths are separate we need 2 BFS’s. 


Additionally, we can use vertices in subgraph K in our longer path besides S and T to 
help us construct a longer path (We will see we can use 2, X and Y, or we can use less, 
but we will never need more). This is because we can make the path from S to X, 
then X to F, then F to Y, then Y to T. So the longer path becomes [S,X,F,Y,T].
The reason we have X and Y is because they are vertexes in subgraph K that are visited by the 2 BFS's 
from F. X can be S. Y can be T. but they 
usally aren't because there are lots of places for the BFS's to hit subgraph K. 


X and Y are needed to make our longer path algorithm polynomial time 
(The end goal is still to find paths [S,F,T] but we want to do this in a smart way that is efficient). 


Naive version of algorithm: 
We will run 2 BFS’s from coordination point F, 
and stop when we visit a vertex in K. 
The first BFS intersects K and identifies a path 
from F to X. The second BFS intersects K and 
identifies a path from F to Y. The longer path is [S,X,F,Y,T]. 
The problem with this is that the graph is directed 
and we cannot use the path from F to X 
(but there may be a path from X to F we can find on the reverse graph, 
so in version 2 we will bfs the reverse graph for this path). 

Version 2:
We will run BFS on the reverse of graph G from 
F, and stop when we intersect K to find the first 
intersection location X. Then we run BFS on normal 
G from coordination point F, and stop when we intersect
K to find Y. The longer path is [S,X,F,Y,T]. This works 
for directed graphs!  However, the problem with this 
version is that there is a possibility that X == Y 
(X cannot equal Y otherwise the longer path is not a simple path). 
Also the X and Y can cause a failure scenario when you attempt 
to create [S,X,F,Y,T] because Y might come topologically 
before X in subgraph K, and in some cases, when Y comes before X, 
the longer path becomes non-simple. WE NEED A NEW VERSION. so lets do version 3. 

Version 3 :
We maintain Xarray and Yarray to find all the locations where coordination point F 
intersects subgraph K. All possible X’s are stored, and all possible Y	‘s 
are stored that can form possible longer paths [S,X,F,Y,T] using F. X’s 
and Y’s are found using the reverse BFS and the normal BFS from version 2. 
We only need one reverse BFS to find all the X’s and one normal BFS to find all the Y’s.  
Then we use a double for loop, to try every possible pair (X,Y) from Xarray 
and Yarray to create longer path [S,X,F,Y,T]. We need to go through every 
possible pair in case of failure scenarios as described in version 2. We 
stop when we find the first valid simple path. This path is a longer path. We are done.
If we don't find a longer path, test the next coordination point F in G-K. 


To check if we are in a failure scenario 
we use merge_two_overlapping_paths_in_dag 
which was also used in LOST EDGES METHOD to merge the paths S->X and Y->T 
in the shortest paths DAG.
Paths X->F and F->Y can also not intersect (they came from seperate BFS's) 
so we make sure those are simple by
using a method similar to merge_two_overlapping_paths_in_dag called:
create_crazy_path_without_overlaps(X, 
                                       Y, 
                                       Z, 
                                       X_to_Z_bfs_tree_parents, 
                                       Z_to_Y_bfs_tree_parents, 
                                       graph, 
                                       shortest_paths_dag_vertices, 
                                       DEBUG=DO_DEBUG) 
which is polynomial.

If we cannot find a longer path from lost edges method, or outer vertex method, 
a longer path does not exist with any set of vertices in G so return False.
There are a lot of optimizations in the original paper. However, 
the following described is a polynomial time algorithm. 

QUICK AND DIRTY RUNTIME ANALYSIS: (we talk about outer vertex method because it is much slower
than lost edges method and dominates the worst case runtime):

We run BFS a lot but never in an exponential way. 
2 BFS’s per vertex and a double for loop for checking every possible pair (X,Y) to 
construct possible paths [S,X,F,Y,T] (checking every possible pair in 
the very worst case is upper bounded by something a lot smaller than V^2). 
Checking failure conditions for one possible [S,X,F,Y,T] is done with 2 DFS’s and 2 BFS's
(done by merge_two_overlapping_paths_in_dag) and 2 DFS's (done by create_crazy_path_without_overlaps).

So worst case for this algorithm is O( 2*(V+E)*V + (6*(V+E) * V^2) * V ) 
which is O((V+E) * V^3), and this is polynomial time.

