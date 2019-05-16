CS 341 Challenge #10: Paths in graphs

Completed by Harman Singh and Vaibhav Khaitan


Consider the following problem: 
given an unweighted directed graph G, and nodes s and t, 
decide whether there exist at least two directed SIMPLE paths from s to t of 
different lengths. Is the problem solvable in polynomial time, or is it NP-complete?

##########################################################################################################################################################

The problem is solvable in polynomial time.
For a brief intuitive understanding of the algorith, read the following: https://docs.google.com/document/d/1TRAFvcLphKYWeRWB8b0ebDvTZB0Q6r1xdhyB3acIMdU/edit?usp=sharing

The following is the algorithm :

We will make a reversed version of the graph G and call this G' 

The algorithm makes use of a subgraph of G we call shortestPathsDAG. 
shortestPathsDAG is a DAG which contains all the shortest paths of equal length from S to T in G.

The number of shortest paths between two nodes in a graph is super exponential.
We can build a structure that holds these paths cheaply, 
avoiding the super exponential runtime. 
shortestPathsDAG will be this structure and we will use it to 
generate any shortest path we desire in linear time.
shortestPathsDAG will start with S, its root, and have multiple paths come out of S, which all end at T. 
Utilizing shortestPathsDAG effectively will allow us to solve the algorithmic challenge 
in polynomial time. Finally, our subgraph called shortestPathsDAG can be built in linear time.
The linear time algorithm is covered below; a quick rundown of the algorithm to build shortestPathsDAG would be the following:
Run BFS on S to get the shortest path from S to every vertex. Run BFS on the reversed graph from T, to get the shortest
path from every vertex to T. Then enumerate all the vertices X in the graph besides S and T, and check if  [S->X->T] is 
the same length as the shortest length path between S and T and if it is, add all of [S->X->T]'s vertices, and edges 
to the shortestPathsDAG. Once all the vertices have been enumerated, shortestPathsDAG contains all the shortest paths between S and T.  

There are 2 ways to create longer paths using shortestPathsDAG:
1) First way is called Lost Edge Method.
   One way is to check every vertex in shortestPathsDAG to see if there is an edge, E, 
   between that vertex and another vertex in shortestPathsDAG 
   that hasn't been added to shortestPathsDAG. (In other words E is not 
   in the set of edges of shortestPathsDAG)
   This edge was not added because it was not in any of the shortest paths 
   we found in the algorithm above to construct shortestPathsDAG. 
   Finding such an edge, we call them LOST_EDGES_IN_THE_DAG or LOST_EDGES, will trivally make a longer path:
   Assume the 2 vertices with the LOST_EDGE is X and Y in the shortestPathsDAG, and the LOST_EDGE goes from X to Y.
   Then the longer path is [S->X, LOST_EDGE, Y->T]. 
   LOST_EDGE was not included in the DAG because it does not help a path 
   get closer from S to T (so cannot be a forward edge in shortestPathsDAG).
   Therefore, LOST_EDGE can only make the path longer, because it is either a cross edge 
   or a backward edge in the shortestPathsDAG. Therefore [S->X, LOST_EDGE, Y->T] is a longer path.
   (The exact algorithm to find a lost edge and use it is 
   discussed in PseudoCode below, and that algorithm is linear time)

2)  Second way is called Outer Vertex Method. 
    A longer path requires vertices to be used outside of shortestPathsDAG.
    To create the longer path we have one very strong condition: 
    It must leave a vertex of the shortestPathsDAG, go through some edges not in the 
    shortestPathsDAG (these edges will be called the CRAZY SEGMENT), 
    and then come back to a vertex in the shortestPathsDAG. The path from S to the ExitPoint (called X), 
    X to the Rentry point back into the shortestPathsDAG (called Y), and Y to T will be a possible long path to consider. 
    There is only 1 Exit Point (we will call this X from now on), 
    and 1 Reentry point(we will case this Y from now on) when we create a longer path 
    (So every long path we consider will have only 1 crazy segment).
    In the proof of correctness we explain why you only need 1 X and 1 Y (and not 2, 3, etc).
    
    We can only find paths of longer length with this approach
    because shortestPathsDAG already contains all the shortest paths between S and T.

How do we discover crazy segments and make longer paths?
All we need to do is find X and Y on the shortestPathsDAG (X and Y are vertices
in the shortestPathsDAG. X and Y are also the intersection points between shortestPathsDAG and the CRAZY SEGMENT).
To do this, we test all the vertices in {G vertices}-{shortestPathsDAG vertices} 
and on each vertex Z, run BFS starting on Z, on the graph(G) to find Y 
(call this CrazyBFS. Store all the Y's you find in an array called possibleYFromZ), 
and run BFS starting on Z, on the reverse of the graph(G') to find X 
(call this ReverseCrazyBFS. Store all the X's you find in an array called possibleXToZ). 

The crazy segment is a path running from one of the X's to Z to one of the Y's. 
We do BFS on the reverse graph to get X, 
because this is the EXIT POINT out of the 
DAG from S to reach Z (The EXIT POINT can also be S in which case S == X). 
We do normal BFS on the graph to get Y, because this will 
find an entry point into the DAG from Z to finally reach T
(The RENTRY POINT can also be T in which case Y == T). 

So the CRAZY SEGMENT is [X -> Z -> Y], and the longer path we can create will be [S->X->Z->Y->T]
The Z we use is called the COORDINATION POINT, because it is where X and Y sprout from in the algorithm. 

(We will later see we wont have to run CrazyBFS and ReverseCrazyBFS 
on every vertex in {G vertices}-{shortestPathsDAG vertices}, 
but actually just a very few number) 

Lets discuss what the {G vertices}-{shortestPathsDAG vertices} 
looks like, call this set subtraction verticesToTest.
verticesToTest are vertices that are in the graph that 
weren't used to construct a path from s to t. We have to use these
vertices, and the edges between them to get a longer path. 
These edges form a forest of components sprawled around the 
shortestPathsDAG. We have to test these components, 
and throw them out if they don't identify a path of longer length.

Finally lets discuss how we create the path given we find an X and Y for a Coordination point Z. 
The path we construct will always have the form [S->X, X->Z, Z->Y, Y->T] 
and X and Y are points in the shortestPathDAG.
Path Creation will always work when X is topologically before Y in the shortestPathDAG. 
If X is topologically after Y, then, for some cases, PATH CREATION will fail.
The is because you would go from S to X, take a path backwards to Y (using crazy segment), and then hope
there is another path in the DAG which can be traversed to get from Y to T, without bumping into the S->X vertices 
in the shortestPathsDAG. In other words, this case will only succeed
when the shortestPathsDAG has multiple branching paths you can traverse to get from S to T!

But now you can understand why we need an array of X's and an array of Y's (possibleXToZ and possibleYFromZ, respectively)
and not just one X and one Y for each Z. We have to try every pair (X,Y) in our Path Creation algorithm 
so that if pair (X,Y) fails we can try another pair and exhaust the possibilities for Z. 
Don't worry our algorithm will still be polynomial time, trust us! Because we will be smart about this!
One technique we will use is memoization for testing these pairs so we aren't retesting
the same X and Y, when we run Outer Vertex Method on a different Z. We will memoize results in 
a map called BADXY = {}, which will map for every vertex X (Exit Vertex), 
the set of Y's that fail for that vertex (Bad Rentry points).
You will see that BADXY gets built very quickly with just one pair (X,Y) 
because the path creation algorithm will do a lot of computations for just one pair.
(The path creation algorithm will do 4 DFS's back to back in the worst case.
 The visited vertices in those 4 DFS's will fill up the BADXY). 

The path creation algorithm is discussed in its helper function in the pseudocode. 

OK LETS PUT ALL THESE PIECES TOGETHER FOR OUTER VERTEX METHOD:

The core idea of pseudo code below: 

Steps: 

while verticesToTest is not Empty:
    Z := verticesToTest.pop() 

    1) Run on Z reverse crazy BFS to find possibleXToZ. Also save visitedXtoZ which is visited vertices from this BFS.
    2) Run on Z crazy BFS to find possibleYFromZ. Also save visitedYFromZ which is visited vertices from this BFS.

    3) If we dont find any X, dont find any Y, or dont find either, then we can reduce verticesToTest getting us closer to algorithm termination: 
       If we dont find any X, then there are no exit points that lead out from the DAG to Z.  
       This also means there are no exit points that lead out from the DAG to the vertices visited by the reverse crazy BFS on Z. 
       (Otherwise, we would have found an X that reaches Z).
        We can remove all the vertices we traversed in the reverse crazy BFS from verticesToTest.  

        Similarly, if we dont find any Y, then there are no re-entry points that lead from Z back into the DAG.
        This also means there are no re-entry points that lead into the DAG using Z's visited vertices.  
        (Otherwise, we would have found a Z that reaches Y).
        We can remove all the vertices we traversed in the crazy BFS from verticesToTest. 

    4) If we found one or more X and one or more Y(So len(possibleXToZ) > 0 and len(possibleYFromZ) > 0):

        4a) For each X in possbleXToZ:
            For each Y in possibleYFromZ:
                if( X != Y and Y is not in BADXY[X]): 
                    The BFSTree touches the shortestPathsDAG in two places, called X and Y in the shortestPathsDAG. 
                    We can CREATE a new path (PLEASE LOOK AT ALGO BELOW TO SEE LONGER PATH CREATION ALGORITHM). 
                    path = createLongerPath(S, T, X, Y, CRAZYSEGMENT, shortestPathsDAG, BADXY)

                    The NEW path is [S->X, X->Z->Y, Y->T], called Q.

                    A) Q's length is longer than the shortestPathLength because all the shortest paths are already in shortestPATHDAG. 
                       We are done. We found 2 simple paths of different length, the shortest one, and [S->X, X->Z->Y, Y->T].
                       RETURN TRUE AND THE TWO PATHS.

                    B) PATH [S->X, X->Z->Y, Y->T] could not be created. 
                       The createLongerPath algorithm memoized the failure results. 
        
        4b) We went through all the pairs (X, Y) and they did not work for Z. They will also not work for vertices in 
            visitedXtoZ or for vertices in visitedYFromZ because these vertices will discover the 
            same or smaller subsets of possibleXToZ and possibleYFromZ. (They cannot discover more, otherwise Z would have discovered it too)
            
            verticesToTest -= visitedXtoZ
            verticesToTest -= visitedYFromZ

       

    If we run out of verticesToTest, we have exhausted all possible longer paths that could have existed, and return False.


######################################################################
Running Time. (Assuming set look up, map lookup, set add is all O(1))

The algorithm does buildShortestPathsDAG at the beginning 
which is 2 BFS calls. and a for loop through every vertex in G in the worst case.
The algorithm also builds shortest paths and puts them in the DAG.
It doesn't iterate through vertices it puts in the DAG. So O( V+E ) runtime. 

The algorithm then checks for lost edges which is O(E*V) time because worst 
case it goes through every vertex and checks every neighbour.

Reverse G which is O(V+E) time. 


In the very worst case, the algorithm does CRAZY BFS AND REVERSE CRAZY BFS 
on every vertex. Each run, verticesToTest goes down by 1.


We go through every pair (X, Y) we find in our 2 BFS's. 
We run the path creation algorithm which is 4 DFS's back to back in the worst case.
So going through each pair in the worst case costs O(V + E).

If we dont consider the heavy optimizations that come from memoization,
then for every vertex, we run through V^2 pairs. 
(if the number of X's we found is V, and number of Y's we found is V)
That means that O( (V+E) * V^3) is the worst case to do outer vertex method. 

I think, with some work, it could be proved that worst case for 
outer vertex method is O(V+E) if we
considered the reduction of verticesToTest 
and memoization of BADXY more thoroughly.  


So worst case for this algorithm is O( 2(V+E)*V + (4(V+E) * V^2) * V ).
which is O((V+E) * V^3)


And this run time is polynomial.


###################################################################################

PSEUDO_CODE:

// THIS IS HELPER FUNCTION 1
// RUNTIME: 2 BFS calls. 
// 
def buildShortestPathsDAG(G, S, T):
    -> Run BFS and get shortest paths from s to every vertex in G. 
        ->Store parentsS[], and distS[] for this BFS.
    -> Reverse the graph and run BFS to get the shortest path from every vertex to T. 
        ->Store parentsT[], and distT[] for this BFS.
    -> Ok store length of shortest path from S to T, call it shortestLength 
        ->(get from one of the BFS's).

       verticesTestedToBeInShortestDAG := All the vertices in the graph except S and T:  
       
       shortestPathsDAG = {} //Adjacency List using Map

       For each vertex X in verticesTestedToBeInShortestDAG: 
            Find distance from vertex X to S (using  DistS[X]). 
            Find distance from vertex X to T (using DistT[X]). 
            Add them to find the distance of path [S->X->T] called length
                if length == shortestLength, 
                    -> construct path [S->X->T] using parentsS and parentsT
                    -> add the path [S->X->T] to the shortestPathsDAG  
                    -> remove all the vertices we added in the path from verticesTestedToBeInShortestDAG

        return {
            "shortestPathsDAG": shortestPatsDAG,
            "shortestPathLength": shortestLength
        }

// THIS IS HELPER FUNCTION 2
// A LOST EDGE is an edge between 2 vertices in the shortestPathDAG that wasn't added to shortestPathDAG
// RUNTIME: Go through every vertex, and check every neighbour worst case.
// So worst case runtime is O(sum of the degree of every vertex in G)  -> O(|E||V|)
def createLongerPathUsingLostEdges(G, shortestPathsDAG):
    
    setOfVerticesInDAG = shortestPathDAG vertices
    // Set look up is fast.

    for each vertex v in setOfVerticesInDAG:
        all_neighbours = G[v]
        dag_neighbours = shortestPathsDAG[v]

        other_neighbours = all_neighbours - dag_neighbours

        for each k in other_neighbours:
            if(k is in setOfVerticesInDAG):
                return True // We found a lost edge. REPORT THAT 2 PATHS OF DIFFERENT SIZE EXIST.
                            // Lost edge creates [S->V, V->K, K->T] which is longer than shortest path 
                            // (check correctness proof for this)
    
    return False //Could not find any lost edges.

//THIS IS HELPER FUNCTION 3. WE USE THIS TO DO CRAZYBFS AND REVERSECRAZYBFS
//The function returns the following: 
// dist => distances from Z to other vertices
// pi => parents 
// visited => the vertices we visited in our BFS
// intersectionVertices => All the vertices we touched in shortestPathsDAG from Z after BFS. 
//                         We only touch the surface of shortestPathDAG and then 
                           break and dont BFS "inside of" shortestPathsDAG
                           If intersection vertices is empty, we did not find an 
                           intersection with DAG.
def crazyBFS(Graph, Z, shortestPathsDAG):  // HELPER FUNCTION DEFINED HERE.
    dist[s] := 0;
    Q := empty queue;
    visited = set()
     
    intersectionVertices = set()
    for all u in V do
        dist[u] := infinity;
        pi[u] := nil;
    
    enqueue(Q, Z);
    while (Q is nonempty) do
        u := dequeue(Q);
        if(u in visited):
            continue
        else: 
            visited.add(u)
        
        if(u in shortestPathsDAG vertices): //This is map look up O(1)
            intersectionVertices.add(u)
            continue //Just keep doing BFS with other vertices. DONT GO INSIDE OF THE DAG after we touch surface.
        
        for all edges (u,v) in Graph do //Out neighbours traversed (Directed BFS)
            enqueue(Q,v)
            dist[v] := dist[u]+1;
            pi[v] := u;


    return ({"dist": dist,
             "pi": pi, 
             "visited": visited, 
             "intersectionVertices": intersectionVertices
             })

// THIS IS HELPER FUNCTION 4. WE USE THIS TO CREATE THE LONGER PATH AFTER FINDING X AND Y.
// X is the exit point of the DAG
// Y is the re-entry point into the DAG
// Z is the coordination point (the point we run crazy bfs and reverse crazy bfs on)
// s is the root of the DAG 
// t is the end of the DAG.
def createLongerPath(s, t, X, Y, CRAZYSEGMENT, shortestPathsDAG, BADXY):

    1)  -> To create the longer path we run DFS to get a path from S->X and we run DFS to get a path from Y->T on the DAG
            -> If Y comes topologically after X, there wont be a problem because CRAZY SEGMENT will be forward edges from X to Y and the DAG 
               points forward from S to T
            -> If Y comes topologically before X, there may be a problem because CRAZY SEGMENT will be backward edges from X to Y 
               so the path from Y->T may intersect with the path from S->X (we cannot have this happen because it must be a SIMPLE Path)
            -> The create path algorithm will assume the Y before X scenerio because it works on the X before Y scenerio:

               1) DFS from S->X on shortestPathsDAG and get a path, and store visited vertices in visitedX. if S->X DFS see's vertex Y, 
                   backtrack and try different path.
                   1a) If DFS fails, and we can't avoid Y, then return False. MEMOIZE FAILURE:
                          for all X in visitedX:
                               BADXY[X].add(Y).

                    1b) We found a path. Go to Step 2.

                2) Then DFS from Y->T on shortestPathsDAG without visiting the vertices from from visitedX or X. 
                    2a) If T is reached, we are done. Return True. 
                    2b) Otherwise T was not reached. Go to Step 3.

                3) Lift restriction that Y->T path can't use visitedX, and try to use visitedX vertices to get to T (this DFS however, still cant visit X). 
                   Save visited vertices to visitedY, including the ones we see when we lift the restriction to not use visitedX. 
                   3a)  If T still cant be reached because Y->T path still sees X, then FAIL PATH CREATION. Return false. MEMOIZE FAILURE:

                        for all X in visitedX:
                            for all Y in visitedY:
                                if(X != Y):
                                    BADXY[X].add(Y).

                   3b)  If T is reached. Go To Step 4. 

                4) We need to do DFS 1 final time on S->X and attempt to get path, without seeing Y or visitedY. Store the dfs vertices in visitedX2.
                           4a) S->X worked. Then Done. Return True.
                           4b) S->X did not work. Return false. The piece of the path we gave up to Y->T was the only segment we could use to go 
                               from S->X without seeing Y or its vertices. Y->T needed this segment too, otherwise, 
                               the previous DFS would have suggested another way for Y to reach T.   
                               Fail Path Creation due to contention for this piece of critical segment in the DAG that both S->X and Y->T needed.
                               MEMOIZE FAILURE:

                               for all X in visitedX2:
                                   for all Y in visitedY:
                                   if(X != Y):
                                        BADXY[X].add(Y).

        -> The above path creation algorithm uses 4 DFS's in the very worst case back to back.
        -> The 4 DFS's can be optimized with early failure stopping by using BADXY (We will talk about that in another paper).
        -> If the memoization above seems confusing to you, we dont need to do it to achieve polynomial time for this problem. 
           just memoize (X,Y) pair instead.



//THIS IS THE MAIN FUNCTION THAT SOLVES THE ALGORITHMIC CHALLENGE:
def two_paths(G, s, t):

    shortestPathDAGBuildResult  = buildShortestPathsDAG(G, s, t) //Call helper function 1
    shortestPathsDAG = shortestPathDAGBuildResult["shortestPathsDAG"]
    shortestPathLength =  shortestPathDAGBuildResult["shortestPathLength"]

    resultOfLostEdgesMethod = createLongerPathUsingLostEdges(G, shortestPathsDAG):
    
    if(resultOfLostEdgesMethod == TRUE):
        return TRUE  //WE ARE DONE, FOUND 2 SIMPLE PATHS OF DIFF LENGTH

    //OTHERWISE DO OUTER VERTEX METHOD.

    ReverseG := Reverse the graph G
    verticesToTest = set() = G vertices - shortestPathsDAG vertices
    BADXY = {} // This maps a vertex X to a set of bad vertices Y that fail path creation 

    while(verticesToTest not empty): // Go through all the vertices in verticesToTest and find X and Y using ReverseCrazyBFS and CrazyBFS
        Z = verticesToTest.pop()

        XResult = crazyBFS(ReverseG, Z, shortestPathsDAG) //We have to BFS on reverse graph to find exit point out of DAG
        YResult = crazyBFS(G, Z, shortestPathsDAG) //We BFS on normal graph to find entry point into DAG
        
        if(len(XResult.intersectionVertices) > 0 and len(YResult.intersectionVertices) > 0): 

           For X in XResult.intersectionVertices:
                For Y in YResult.intersectionVertices:
                     if( X != Y and Y is not in BADXY[X]):
                        // PATH IS [S->X->Z->Y->T]
                        CRAZYSEGMENT = Create Path From X to Z to Y 
                                         -> Create using ancestors of X stored in XResult.pi to Z.
                                         -> Create using ancestors of Y stored in YResult.pi to Z.
                                         -> Combine 2 and that is CrazySegment.

                        NEWPATH_RESULT = createLongerPath(s, t, X, Y, CRAZYSEGMENT, shortestPathsDAG, BADXY)
                        
                        if NEWPATH_RESULT == TRUE:
                            RETURN TRUE (2 PATHS ARE [S->X->Z->Y->T] AND shortest STPATH)
            // All pairs failed
            verticesToTest -= XResult.visited
            verticesToTest -= YResult.visited;
                            
        else: //Either no X could be found, or no Y could be found, or neither could be found. 
            if(len(XResult.intersectionVertices) == 0):
                // X could not be found which means no exit point out of DAG 
                // All the vertices we saw in the BFS will not be able to find exit points either
                // so remove all of them as potential COORDINATION POINTS
                verticesToTest -= XResult.visited

            if(len(YResult.intersectionVertices) == 0):
                // Y could not be found which means no REENTRY point into the DAG 
                // All the vertices we see also cannot find a REENTRY POINT INTO DAG 
                // so remove all of them as potential COORDINATION POINTS
                verticesToTest -= YResult.visited; //REMOVE ALL VISITED VERTICES.         

    RETURN FALSE

###################################################################################
