
import numpy as np
import pprint
from collections import defaultdict, deque
from utility import bfs, get_path_to_root, create_example_rand_directed_graph


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
def create_shortest_paths_dag(graph, reversed_graph, s, t):
    # You can do this by doing BFS on every vertex. Then
    
    graphBFS = bfs(graph, s)
    reverseGraphBFS = bfs(reversed_graph, t)
    
    parentsS = graphBFS["parent"]
    distS = graphBFS["dist"]
    parentsT = reverseGraphBFS["parent"]
    distT = reverseGraphBFS["dist"]
    
    shortestLength = distS[t] # also equal to distT[s]

    print("shortest_path_length is ", shortestLength)

    # test all vertices except s and t to be in the subgraph
    vertices_to_test = set(graph.keys()) - set([s, t]) 
    print("dist S", distS)
    print("dist T",  distT)

    shortest_path_dag = defaultdict(set)
    shortest_path_vertices = set() # have to maintain seperately

    print(vertices_to_test)

    for v in vertices_to_test:
        print("v is ", v)
        dist_from_s_to_v = distS.get(v)
        dist_from_v_to_t = distT.get(v)
        
        # If the distances exist because they were traversed then
        if(dist_from_s_to_v and dist_from_v_to_t):
            length = dist_from_s_to_v + dist_from_v_to_t # Add them to find distance of path [S -> v -> T]
            print("length was",  length)
            
            if(length == shortestLength):
                # Add path [s -> v -> t] to shortest paths subgraph
                path_V_to_S = get_path_to_root(parentsS, v)
                path_V_to_T = get_path_to_root(parentsT, v) # remove the V node
                print("path V to S", path_V_to_S)
                print("path V to T", path_V_to_T)

                path_S_to_V_to_T = path_V_to_S[::-1] + path_V_to_T[1::]
                
                print("path S to V to T", path_S_to_V_to_T)
                # add path to subgraph
                 
                add_path_to_graph(path_S_to_V_to_T, shortest_path_dag)
                shortest_path_vertices.update(path_S_to_V_to_T)

                # TODO: you can optimize here by subtracting vertices you added from vertices to test. 
                             
    print("shortest path subgraph is: ", shortest_path_dag)
    print("graph_bfs_tree_with_root_s", graphBFS)
    print("reverse_graph_bfs_tree_with_root_t", reverseGraphBFS)
    return {
        "shortest_path_dag": shortest_path_dag,
        "shortest_path_dag_vertices": shortest_path_vertices, 
        "a_shortest_path": get_path_to_root(parentsS, t)[::-1],
        "shortest_path_length": shortestLength,
        "graph_bfs_tree_with_root_s": graphBFS["parent"],
        "reverse_graph_bfs_tree_with_root_t": reverseGraphBFS["parent"]
    }



def crazy_bfs(graph, Z, shortest_paths_dag, vertices_in_dag):
    seen, queue = set([Z]), deque([Z])
    parent, dist = {}, {}
    dist[Z] = 0
    parent[Z] = None # root has no parent

    intersection_vertices = set()

    while queue:
        v = queue.popleft()

        if(v in seen):
            continue
        else:
             seen.add(node)

        if(v in vertices_in_dag):
            intersection_vertices.add(v)
            continue
            # KEEP DOING BFS WITH OTHER VERTICES. DONT GO INSIDE THE DAG after we touch surface

        for node in graph[v]: 
            dist[node] = dist[v] + 1
            parent[node] = v
            queue.append(node)
    
    return ({
        "seen": seen,
        "parent": parent,
        "dist": dist,
        "intersection_vertices": intersection_vertices,
    })


def create_longer_path_using_lost_edges(graph, 
                                        shortest_path_dag, 
                                        vertices_in_dag, 
                                        graph_bfs_tree_with_root_s, 
                                        reverse_graph_bfs_tree_with_root_t):
    for V in vertices_in_dag: 
        all_neighbors = graph[V]
        dag_neighbors = shortest_path_dag[V]

        not_in_dag_neighbors = all_neighbors - dag_neighbors # set subtraction
        
        for K in not_in_dag_neighbors:
            if(K in vertices_in_dag):
                # We found a lost edge. REPORT THAT 2 PATHS OF DIFFERENT SIZE EXIST.
                # Lost edge creates [S->V, V->K, K->T] which is longer than shortest path, 
                # because otherwise that path would be in the shortest_path_dag
                # Question why would K -> T exist? it has to exist because K is a vertice in the DAG
                return {
                    "result": True,
                    "a_longer_path": get_path_to_root(graph_bfs_tree_with_root_s, V)[::-1] + [V, K] + get_path_to_root(reverse_graph_bfs_tree_with_root_t, K)
                    }
    
    return {
        "result": False,
        "a_longer_path": None
        }


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
# WE USE 4 DFS's to exhaust these possibilities in the shortest paths dag
# and ensure that [S->X->Z->Y->T] is a simple path.

def create_longer_path_using_outer_vertex(s, t, X, Y, shortest_paths_dag):
    '''
     PSUEDOCODE FOR THE 4 DFS's in this function numbered 1 to 4:

     1) DFS from S->X on shortestPathsDAG and get a path, and store visited vertices in visitedX. 
        if S->X DFS see's vertex Y,  backtrack and try different path.
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
    pass

def solution(graph, s, t):

    # Reversed graph is needed. 
    reversed_graph = defaultdict(set) # a map with values as list type

    # Go through graph and create reversed graph
    for parent, neighbors in graph.items():
        for n in neighbors:
            reversed_graph[n].add(parent)
    
    pprint.pprint(reversed_graph)

    buildResult = create_shortest_paths_dag(graph, reversed_graph, s, t)

    # VARIABLE BELOW IS THE REASON WHY ALGO WILL BE POLYNOMIAL
    shortest_path_dag = buildResult["shortest_path_dag"] # DIRECTED ACYCLIC GRAPH that contains all shortest paths from s to t,
    shortest_path_dag_vertices = buildResult["shortest_path_dag_vertices"] # VERTICES OF THE DAG ABOVE
    shortest_path_length = buildResult["shortest_path_length"] # shortest path from s to t

    # BELOW VARIABLES NOT USED FOR ALGORITHM. JUST USED TO CREATE THE 2 SIMPLE PATHS ONCE WE FIND THEM!!
    # THEY are used for instance in create_longer_path_using_lost_edges only when 2 simple paths had been found 
    a_shortest_path = buildResult["a_shortest_path"]
    graph_bfs_tree_with_root_s = buildResult["graph_bfs_tree_with_root_s"] 
    reverse_graph_bfs_tree_with_root_t = buildResult["reverse_graph_bfs_tree_with_root_t"]
    
    # WE TRY 2 METHODS to create our 2 paths: LOST_EDGES method, and OUTER VERTEX Method. 
    # if both methods fail, THERE IS NO 2 PATHS, where 1 is longer than the other from S to T 
    # These 2 methods encapsulate all possible ways for these 2 paths to exist, that's why if they fail, then there is none

    did_we_get_2_paths = create_longer_path_using_lost_edges(graph, 
                                                             shortest_path_dag, 
                                                             shortest_path_dag_vertices, 
                                                             graph_bfs_tree_with_root_s, 
                                                             reverse_graph_bfs_tree_with_root_t)

    if(did_we_get_2_paths["result"] == True):
        # WE FOUND 2 PATHS!!!
        print("WE FOUND A SHORTER AND A LONGER DIRECTED SIMPLE PATH USING LOST EDGES METHOD. THEY ARE: ")
        print(a_shortest_path)
        print(did_we_get_2_paths["a_longer_path"])

    # OK LOST EDGES FAILED. now have to try outer vertex method.

    vertices_to_test = set(g.keys()) - shortest_path_dag_vertices # possible coordinate vertices
    
    # Crazy bfs is done to find places where we touch shortest_path_dag from vertices not in the dag
    # this is why its called the outer vertex method. we use a vertex on the outside to go into the shortest_paths_dag
    # we are attempting to make a longer st-path by doing the following: [S->X->Z->Y->T]. 
    # X, and Y are places where we touch the shortest_paths_dag from coordinate vertex Z (Z is from vertices_to_test variable above)
    while vertices_to_test:
        Z = vertices_to_test.pop() # Get a random vertex from set and do crazy bfs and reverse crazy bfs on it.
        
        # We have to BFS on reverse graph to find exit point out of DAG
        XResult = crazy_bfs(reversed_graph, Z, shortest_path_dag, shortest_path_dag_vertices)
        # We BFS on normal graph to find entry point into DAG
        YResult = crazy_bfs(graph, Z, shortest_path_dag, shortest_path_dag_vertices)

        # WE CAN ADD MORE DYNAMIC PROGRAMMING HERE TO MAKE THIS VERY FAST. 
        # but we wont for now because without it, its still polynomial time. check readme for extra dp.

        for x in XResult["intersection_vertices"]:
            for y in YResult["intersection_vertices"]:
                if( x != y ): 
                    # Possible path can be [S->X->Z->Y->T]
                    # However we must check for some bad cases
                    pass



###############################################################################################
############################################################################################

# EXECUTE POLYNOMIAL TIME SOLUTION TO PROBLEM HERE:
g = create_example_rand_directed_graph(10)

pprint.pprint(g)

# DOES THERE EXIST AT LEAST 2 DIRECTED SIMPLE PATHS FROM S TO T of different lengths
solution(graph=g, s=1, t=5)
    



