
import numpy as np
import pprint
from collections import defaultdict, deque
from utility import bfs, dfs_with_restriction_set, get_path_to_root, create_example_rand_directed_graph


'''
Consider the following problem: 
given an unweighted directed graph G, and nodes s and t, 
decide whether there exist at least two directed SIMPLE paths from s to t of 
different lengths. Is the problem solvable in polynomial time, or is it NP-complete?


THE PROBLEM IS SOLVABLE IN POLYNOMIAL TIME. SOLUTION BELOW. DONE WITH LOTS OF BFS, DFS, AND USING A DIRECTED 
ACYCLIC GRAPH AS A MAP STRUCTURE TO MEMOIZE LOTS OF SOLUTIONS QUICKLY TO ELIMINATE THE SUPER-EXPONENTIAL TIME. 

KEY LEARNING: DIRECTECD ACYCLIC GRAPHS ARE USEFUL FOR DYNAMIC PROGRAMMING when the problem has to do with enumerating 
all possible paths in graphs. 
'''


DO_DEBUG=True


##### START READING THE CODE FROM HERE. THIS IS THE MAIN METHOD. 
def poly_solution(graph, s, t, DEBUG=DO_DEBUG):

    # Reversed graph is needed. 
    reversed_graph = defaultdict(set) # a map with values as list type

    # Go through graph and create reversed graph
    for parent, neighbors in graph.items():
        for n in neighbors:
            reversed_graph[n].add(parent)
    
    pprint.pprint(reversed_graph)

    buildResult = create_shortest_paths_dag(graph, reversed_graph, s, t)

    if(buildResult["is_there_shortest_path"] == False):
        if(DEBUG): print("THERE IS NO SHORTEST PATH BETWEEN S AND T, SO THERE IS NO SOLUTION. FALSE. BYE")
        return {
            "result": False,
           
        }


    # VARIABLE BELOW IS THE REASON WHY ALGO WILL BE POLYNOMIAL
    shortest_paths_dag = buildResult["shortest_paths_dag"] # DIRECTED ACYCLIC GRAPH that contains all shortest paths from s to t,
    shortest_path_length = buildResult["shortest_path_length"] # shortest path from s to t

    # BELOW VARIABLES NOT USED FOR ALGORITHM. JUST USED TO CREATE THE 2 SIMPLE PATHS ONCE WE FIND THEM!!
    # THEY are used for instance in create_longer_path_using_lost_edges only when 2 simple paths had been found 
    a_shortest_path = buildResult["a_shortest_path"]

    if(DEBUG): print("A SHORTEST PATH IS: ", a_shortest_path)
    if(DEBUG): print("SHORTEST PATH LENGTH IS ", shortest_path_length)

    graph_bfs_tree_with_root_s = buildResult["graph_bfs_tree_with_root_s"] 
    reverse_graph_bfs_tree_with_root_t = buildResult["reverse_graph_bfs_tree_with_root_t"]
    
    # WE TRY 2 METHODS to create our 2 paths: LOST_EDGES method, and OUTER VERTEX Method. 
    # if both methods fail, THERE IS NO 2 PATHS, where 1 is longer than the other from S to T 
    # These 2 methods encapsulate all possible ways for these 2 paths to exist, that's why if they fail, then there is none

    did_we_get_2_paths_using_lost_edges = create_longer_path_using_lost_edges(graph, 
                                                             s,
                                                             t,
                                                             shortest_paths_dag, 
                                                             graph_bfs_tree_with_root_s, 
                                                             reverse_graph_bfs_tree_with_root_t)

    if(did_we_get_2_paths_using_lost_edges["result"]):
        # WE FOUND 2 PATHS!!!
        if(DEBUG): print("WE FOUND A SHORTER AND A LONGER DIRECTED SIMPLE PATH USING LOST EDGES METHOD. THEY ARE: ")
        if(DEBUG): print(a_shortest_path)
        if(DEBUG): print(did_we_get_2_paths_using_lost_edges["a_longer_path"])
        
        return {
            "result": True,
            "a_shortest_path": a_shortest_path,
            "a_longer_path": did_we_get_2_paths_using_lost_edges["a_longer_path"]
        }

    # OK LOST EDGES FAILED. now have to try outer vertex method.
    did_we_get_2_paths_using_outer_an_outer_vertex = create_longer_path_using_an_outer_vertex(graph=graph, 
                                                                                              reversed_graph=reversed_graph, 
                                                                                              shortest_paths_dag=shortest_paths_dag, 
                                                                                              s=s, 
                                                                                              t=t)   


    if(did_we_get_2_paths_using_outer_an_outer_vertex["result"]):
        if(DEBUG): print("WE FOUND A SHORTER AND A LONGER DIRECTED SIMPLE PATH USING OUTER VERTEX METHOD. THEY ARE: ")
        if(DEBUG): print(a_shortest_path)
        if(DEBUG): print(did_we_get_2_paths_using_outer_an_outer_vertex["a_longer_path"])
        
        return {
            "result": True,
            "a_shortest_path": a_shortest_path,
            "a_longer_path": did_we_get_2_paths_using_outer_an_outer_vertex["a_longer_path"]
        }

    return {
        "result": False
    }



# Add a path, for example path [a,b,c], to a graph g
# graph will then contain a->b->c
def add_path_to_graph(path, g):
    
    # add the path to the graph
    parent = path[0]

    for i in range(1, len(path)):
        child = path[i]
        g[parent].add(child)
        parent = child
    
    # check if last child in graph, if not, add with its value being empty set
    if(path[-1] not in g):
        g[path[-1]] = set()


# Shortest paths subgraph contains all the shortest paths from node s to t.
# it is a DAG from s to t. 
def create_shortest_paths_dag(graph, reversed_graph, s, t, DEBUG=DO_DEBUG):
    # You can do this by doing BFS on every vertex. Then
    
    graphBFS = bfs(graph, s)
    reverseGraphBFS = bfs(reversed_graph, t)
    
    parentsS = graphBFS["parents"]
    distS = graphBFS["dist"]
    parentsT = reverseGraphBFS["parents"]
    distT = reverseGraphBFS["dist"]
    
    shortestLength = distS.get(t) # also equal to distT[s]

    if(shortestLength is None):
        # There is no path from s to t, FAIL!
        return {
            "is_there_shortest_path": False
        }

    if(DEBUG): print("shortest_path_length is ", shortestLength)
    if(DEBUG): print("shortest length other way is ", distT[s])

    # test all vertices except s and t to be in the subgraph
    vertices_to_test = set(graph.keys()) - set([s, t]) 
    if(DEBUG): print("dist S", distS)
    if(DEBUG): print("dist T",  distT)

    shortest_paths_dag = defaultdict(set)
    # shortest_path_vertices = set() # have to maintain seperately

    if(DEBUG): print(vertices_to_test)

    for v in vertices_to_test:
        if(DEBUG): print("v is ", v)
        dist_from_s_to_v = distS.get(v)
        dist_from_v_to_t = distT.get(v)
        
        # If the distances exist because they were traversed then
        if(dist_from_s_to_v and dist_from_v_to_t):
            length = dist_from_s_to_v + dist_from_v_to_t # Add them to find distance of path [S -> v -> T]
            if(DEBUG): print("length was",  length)
            
            if(length == shortestLength):
                # Add path [s -> v -> t] to shortest paths subgraph
                path_V_to_S = get_path_to_root(parentsS, v)
                path_V_to_T = get_path_to_root(parentsT, v) # remove the V node
                if(DEBUG): print("path V to S", path_V_to_S)
                if(DEBUG): print("path V to T", path_V_to_T)

                path_S_to_V_to_T = path_V_to_S[::-1] + path_V_to_T[1::]
                
                if(DEBUG): print("path S to V to T", path_S_to_V_to_T)
                # add path to subgraph
                 
                add_path_to_graph(path_S_to_V_to_T, shortest_paths_dag)
                # shortest_path_vertices.update(path_S_to_V_to_T)

                # TODO: you can optimize here by subtracting vertices you added from vertices to test. 
                             
    if(DEBUG): print("shortest path subgraph is: ", shortest_paths_dag)
    if(DEBUG): print("graph_bfs_tree_with_root_s", graphBFS)
    if(DEBUG): print("reverse_graph_bfs_tree_with_root_t", reverseGraphBFS)
    return {
        "is_there_shortest_path": True,
        "shortest_paths_dag": shortest_paths_dag,
        "a_shortest_path": get_path_to_root(parentsS, t)[::-1],
        "shortest_path_length": shortestLength,
        "graph_bfs_tree_with_root_s": graphBFS["parents"],
        "reverse_graph_bfs_tree_with_root_t": reverseGraphBFS["parents"]
    }



def crazy_bfs(graph, Z, shortest_paths_dag, DEBUG=DO_DEBUG):
    seen, queue = set(), deque([Z])
    parent, dist = {}, {}
    dist[Z] = 0
    parent[Z] = None # root has no parent

    vertices_in_dag = set(shortest_paths_dag.keys())
    
    intersection_vertices = set()

    while queue:
        v = queue.popleft()

        if(v in seen):
            continue
        else:
             seen.add(v)

        if(v in vertices_in_dag):
            intersection_vertices.add(v)
            continue
            # KEEP DOING BFS WITH OTHER VERTICES. DONT GO INSIDE THE DAG after we touch surface

        for node in graph[v]: 
            dist[node] = dist[v] + 1
            parent[node] = v
            queue.append(node)
    
    if DEBUG: print("CRAZY BFS GRAPH FOR " + str(Z) + " is the following: " + str(graph))
    if DEBUG: print("CRAZY BFS SEEN FOR " + str(Z) + " is the following: " + str(seen))
    if DEBUG: print("CRAZY BFS PARENTS FOR " + str(Z) + " is the following: " + str(parent))
    if DEBUG: print("CRAZY BFS DIST FOR " + str(Z) + " is the following: " + str(dist))
    if DEBUG: print("CRAZY BFS INTERSECTION VERTICES FOR " + str(Z) + " is the following: " + str(intersection_vertices))
 
    return ({
        "seen": seen,
        "parents": parent,
        "dist": dist,
        "intersection_vertices": intersection_vertices,
    })


# This method merges overlapping paths together in the shortest_paths_dag
# It finds two paths in the DAG, one that is s->X and one that is Y->t, and these two paths DONT OVERLAP!
# Overlapping occurs due to Y topologically coming before X in the DAG so the two paths created overlap
# 
# This method tries to find 2 paths that go from s->X and Y->t without either overlapping to create SIMPLE paths 
# When the function returns, it sends back 2 parent arrays, one for each path, to allow for the creation of the 
# longer path (or fails because every possible 2 paths will overlap when S is trying to reach X and Y is trying to reach T in their respective DFS's

# This functions is used by Lost edges method and outer vertex method which merge the two paths returned either with a 
# lost edge (edge between two vertices in shortest paths dag, but edge isnt in shortest paths dag. edge is in graph however!), 
# or outer vertex (an outer vertex is a vertex in the graph but not in the DAG) 
# respectively to create the longer path from s->t
# RUNTIME ANALYSIS: DFS on a DAG is O(V+E). WE DO 4 IN THIS BACK TO BACK WORST CASE SO 4 * O(V+E)
def merge_two_overlapping_paths_in_dag(s, t, X, Y, shortest_paths_dag, DEBUG=DO_DEBUG):
    '''
     PSUEDOCODE FOR THE 4 DFS's in this function numbered 1 to 4:

     1) DFS from S->X on shortestPathsDAG and get a path, and store visited vertices in visitedX. 
        if S->X DFS see's vertex Y,  give up on path, continue dfs search with other nodes.
        1a) If DFS fails, and we can't avoid Y, then return False. 
        1b) We found a path. Go to Step 2
     
     2) Then DFS from Y->T on shortestPathsDAG without visiting the vertices from from visitedX or X. 
         2a) If T is reached, we are done. Return True. 
         2b) Otherwise T was not reached. Go to Step 3

     3) Lift restriction that Y->T path can't use visitedX, and try to use visitedX vertices to get to T (this DFS however, still cant visit X). 
        Save visited vertices to visitedY, including the ones we see when we lift the restriction to not use visitedX. 
        3a)  If T still cant be reached because Y->T path still sees X, then FAIL PATH CREATION. Return false. 
        3b)  If T is reached. Go To Step 4
     4) We need to do DFS 1 final time on S->X and attempt to get path, without seeing Y or visitedY. Store the dfs vertices in visitedX2.
                4a) S->X worked. Then Done. Return True.
                4b) S->X did not work. Return false. The piece of the path we gave up to Y->T was the only segment we could use to go 
                    from S->X without seeing Y or its vertices. Y->T needed this segment too, otherwise, 
                    the previous DFS would have suggested another way for Y to reach T.   
                    Fail Path Creation due to contention for this piece of critical segment in the DAG that both S->X and Y->T needed.
    '''
    
    S_to_X_dfs = dfs_with_restriction_set(graph=shortest_paths_dag, 
                                          start=s, 
                                          end=X, 
                                          restriction_set=set([Y]))

    if(S_to_X_dfs["result"] == False):
        if(DEBUG): print("FAILED ON FIRST DFS.")
        return {
            "result": False
        }
    
    # Do second dfs

    Y_to_T_dfs = dfs_with_restriction_set(graph=shortest_paths_dag, 
                                          start=Y, 
                                          end=t, 
                                          restriction_set=S_to_X_dfs["seen"])

    if(Y_to_T_dfs["result"] == True):
        if(DEBUG): print("IT WAS MERGED WITHIN 2 DFS'S")
        return {
            "result": True,
            "S_to_X_dfs_tree": S_to_X_dfs["parents"],
            "Y_to_T_dfs_tree": Y_to_T_dfs["parents"] 
        }

    # The second dfs failed, so move on to the third dfs:

    Y_to_T_dfs_2 =  dfs_with_restriction_set(graph=shortest_paths_dag, 
                                          start=Y, 
                                          end=t, 
                                          restriction_set=set([X]))
    if(Y_to_T_dfs_2["result"] == False):
        if(DEBUG): print("FAILED ON THIRD DFS")
        return {
            "result": False
        }
    
    # do fourth dfs

    S_to_X_dfs_2 = dfs_with_restriction_set(graph=shortest_paths_dag, 
                                          start=s, 
                                          end=X, 
                                          restriction_set=Y_to_T_dfs_2["seen"])

    if(S_to_X_dfs_2["result"] == True):
        if(DEBUG): print("IT WAS MERGED WITHIN 4 DFS'S")

        return {
            "result": True,
            "S_to_X_dfs_tree": S_to_X_dfs_2["parents"],
            "Y_to_T_dfs_tree": Y_to_T_dfs_2["parents"] 
        }
    else:
        return {
            "result": False
        }


# We find a backedge in the nodes in shortest path dag, and then try to make a simple path with it. 
# Only uses nodes within the shortest path dag to create a longer path
def create_longer_path_using_lost_edges(graph, 
                                        s,
                                        t,
                                        shortest_paths_dag, 
                                        graph_bfs_tree_with_root_s, 
                                        reverse_graph_bfs_tree_with_root_t, 
                                        DEBUG=DO_DEBUG):
    if(DEBUG): print("Start lost edges method")
    vertices_in_dag = set(shortest_paths_dag.keys())
    for V in vertices_in_dag: 
        all_neighbors = graph[V]
        dag_neighbors = shortest_paths_dag[V]

        # lost dag neighbors are vertices in the dag that have relationships 
        # with other vertices in the dag, but not dag relationships 
        # because dag relationships required a shortest path to exist. 
        # Instead these are relationships that had existed in the original graph
        # since lost dag neighbor relationships didnt make the cut to be part of the shortest_paths_dag, they
        # could be used to make a longer path! But we have to also check that the longer path is a simple path
        # which is done by merge_two_overlapping_paths_in_dag
        lost_dag_neighbors = vertices_in_dag.intersection(all_neighbors - dag_neighbors) # set subtraction
        
        
        if(DEBUG): print("NOT IN DAG NEIGHBORS IS ") 
        if(DEBUG): print(lost_dag_neighbors)
        
        for K in lost_dag_neighbors:            
            # We found a lost edge. REPORT THAT 2 PATHS OF DIFFERENT SIZE EXIST.
            # Lost edge creates [S->V, V->K, K->T] which is longer than shortest path, 
            # V->K is the lost edge. It has length 1 trivially because its an edge. 
            # [S->V, V->K, K->T] is also [S->V] + [K->T] which we can use in merge_two_overlapping_paths_in_dag
            # You have to check if the path is a simple path 
            # How to check simple, why would it not be simple? 
            # Because V->K is a backedge, path S->V and path K->T can overlap
            # one way to fix that is exhuastively checking all S->V, and K->T so that they do not overlap, and 
            # that is a longer path
            # a better way to do this is using merge_two_overlapping_paths_in_dags which only requires 4 DFS's to do this

            did_path_merge_work = merge_two_overlapping_paths_in_dag(s=s, t=t, X=V, Y=K, shortest_paths_dag=shortest_paths_dag)

            if(did_path_merge_work["result"]):
                if(DEBUG): print("S_TO_X_DFS_TREE", did_path_merge_work["S_to_X_dfs_tree"])
                if(DEBUG): print("Y_TO_T_DFS_TREE", did_path_merge_work["Y_to_T_dfs_tree"])

                return {
                    "result": True,
                    "a_longer_path": get_path_to_root(graph_bfs_tree_with_root_s, V)[::-1] + get_path_to_root(reverse_graph_bfs_tree_with_root_t, K)
                    }
    
    if(DEBUG): print("Lost edges method yielded no results")
    return {
        "result": False,
        "a_longer_path": None
        }


    
def create_longer_path_using_an_outer_vertex(graph, reversed_graph, shortest_paths_dag, s, t, DEBUG=DO_DEBUG):
    if DEBUG: print("START OUTER VERTEX METHOD")
    shortest_paths_dag_vertices = set(shortest_paths_dag.keys()) 

    vertices_to_test = set(graph.keys()) - shortest_paths_dag_vertices # possible coordinate vertices
    
    # Crazy bfs is done to find places where we touch shortest_paths_dag from vertices not in the dag
    # this is why its called the outer vertex method. we use a vertex on the outside to go into the shortest_paths_dag
    # we are attempting to make a longer st-path by doing the following: [S->X->Z->Y->T]. 
    # X, and Y are places where we touch the shortest_paths_dag from coordinate vertex Z (Z is from vertices_to_test variable above)
    
    if DEBUG: print("OUTER VERTICES TO TEST", vertices_to_test)

    while vertices_to_test:
        Z = vertices_to_test.pop() # Get a random vertex from set and do crazy bfs and reverse crazy bfs on it.
        
        # We have to BFS on reverse graph to find exit point out of DAG
        X_to_Z_Result = crazy_bfs(reversed_graph, Z, shortest_paths_dag)
        # We BFS on normal graph to find entry point into DAG
        Z_to_Y_Result = crazy_bfs(graph, Z, shortest_paths_dag)

        # WE CAN ADD MORE DYNAMIC PROGRAMMING HERE TO MAKE THIS VERY FAST. 
        # but we wont for now because without it, its still polynomial time. check readme for extra dp.
        if(DEBUG): print("TEST NEW COORDINATION POINT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        if(DEBUG): print("COORDINATION POINT Z IS ", Z)
        if DEBUG: print("CRAZY BFS RESULT X_to_Z_result", X_to_Z_Result)
        if DEBUG: print("CRAZY_BFS_RESULT_Z_to_Y",Z_to_Y_Result)

        for x in X_to_Z_Result["intersection_vertices"]:
            for y in Z_to_Y_Result["intersection_vertices"]:

                if DEBUG: print("WE WILL TEST X, Y PAIR: ", (x,y))

                if( x != y ): 
                    # Possible path can be [S->X->Z->Y->T]
                    # However we must check for some bad cases
                    
                    # Path we want to make is  [S->X->Z->Y->T], 
                    # X is the exit point of the DAG
                    # Y is the re-entry point into the DAG
                    # Z is the coordination point (the point we run crazy bfs and reverse crazy bfs on)
                    # s is the root of the DAG 
                    # t is the end of the DAG.

                    # [S->X->Z->Y->T] can fail if Y comes topologically before X in shortest_paths_dag because 
                    # this will cause a NONSIMPLE PATH (Y->T may intersect with S->X) 
                    # we have to exhaust all possible paths to find segments Y->T and S->X that dont intersect, 
                    # otherwise tell them not possible.  
                    # WE USE 4 DFS's to exhaust these possibilities in merge_two_overlapping_paths_in_dag
                    # and ensure that [S->X->Z->Y->T] is a simple path.
                    
                    longer_path_result = merge_two_overlapping_paths_in_dag(s, t, x, y, shortest_paths_dag)
                    if DEBUG: print("LONGER PATH RESULT FOR (z,x,y) = " + str((Z, x, y)) + "is the following: " + str(longer_path_result))
                    if DEBUG: print("CRAZY BFS RESULT X_to_Z_result[parents]", X_to_Z_Result["parents"])
                    if DEBUG: print("CRAZY_BFS_RESULT_Z_to_Y[parents]",Z_to_Y_Result["parents"])

                    if(longer_path_result["result"]):
                        print("There are is a shorter and longer path! They are the following: ")
                        

                        

                        # create [S->X->Z->Y->T]
                        if(DEBUG): print("S->X", get_path_to_root(longer_path_result["S_to_X_dfs_tree"], x) ) # S is the root of this. get path to s. WE want away from S so flip 
                        if(DEBUG): print("X->Z", get_path_to_root(X_to_Z_Result["parents"], x) ) # Z is the root of this parents array
                        if(DEBUG): print("Z->Y", get_path_to_root(Z_to_Y_Result["parents"], y)[::-1] ) # Z is the root of this parents array
                        if(DEBUG): print("Y->T", get_path_to_root(longer_path_result["Y_to_T_dfs_tree"], t)[::-1]) # Y is the roof of this. Dont want path to Y, but path to T, so flip it with [::-1]

                        a_longer_path =  get_path_to_root(longer_path_result["S_to_X_dfs_tree"], x) + \
                                          get_path_to_root(X_to_Z_Result["parents"], x)[1:] +  \
                                          get_path_to_root(Z_to_Y_Result["parents"], y)[::-1][1:] + \
                                          get_path_to_root(longer_path_result["Y_to_T_dfs_tree"], t)[::-1][1:]

                        if(DEBUG): print("a longest path: "  + str(a_longer_path))
                        return {
                                "result": True,
                                "a_longer_path": a_longer_path
                                }

    if(DEBUG): print("END TEST FOR COORDINATION POINT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


    if DEBUG: print("OUTER VERTEX METHOD DID NOT YIELD RESULTS")

    return {
        "result": False
    }



###############################################################################################
############################################################################################

# EXECUTE POLYNOMIAL TIME SOLUTION TO PROBLEM HERE:
# g = create_example_rand_directed_graph(10)

# pprint.pprint(g)

# DOES THERE EXIST AT LEAST 2 DIRECTED SIMPLE PATHS FROM S TO T of different lengths
# poly_solution(graph=g, s=1, t=5)
    



