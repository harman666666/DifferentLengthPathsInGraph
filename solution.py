
import numpy as np
import pprint
from collections import deque, defaultdict

# Creates random directed unweighted graph. 
def create_example_rand_directed_graph(vertices):
 
    # This is how we do it with 
    # an adjacency matrix representation
    # adjacency = np.random.randint(0,2,(vertices,vertices))

    g = defaultdict(set) # default dict is a key-value map datatype with values that have set type.
    
    # Now create adjacency list representation from this. Map is adjacency list
    for i in range(vertices):
        ne = np.random.randint(vertices) 
        g[i] = set(np.random.choice(vertices, ne, replace=False))


    return g

g = create_example_rand_directed_graph(10)

pprint.pprint(g)


def bfs(graph, s, t):

    seen, queue = set([s]), deque([s])
    parent, dist = {}, {}
    dist[s] = 0
    parent[s] = None
    
    while queue:
        v = queue.popleft()

        for node in graph[v]:
            if(node not in seen):
                seen.add(node)
                dist[node] = dist[v] + 1
                parent[node] = v
                queue.append(node)
    
    return ({
        "seen": seen,
        "parent": parent,
        "dist": dist,
        "shortestPathDistance": dist[t]
    })


result = bfs(g, 1, 3)
pprint.pprint(result)


def get_path_to_root(bfs_tree, node):
    # root has parent called undefiend
    n = node 
    path = [n]

    while True:
        # print("n is", n)
        parent = bfs_tree[n]
        
        if(parent is None):
            break

        path.append(parent)
        n = parent

    return path

def add_path_to_graph(path, g):
    
    # add the path to the graph
    parent = path[0]

    for i in range(1, len(path)):
        child = path[i]
        g[parent].add(child)
        parent = child




# Shortest paths subgraph contains all the shortest paths from node s to t.
# it is a DAG from s to t. 
def create_shortest_paths_subgraph(graph, reversed_graph, s, t):
    # You can do this by doing BFS on every vertex. Then
    
    graphBFS = bfs(graph, s, t)
    reverseGraphBFS = bfs(reversed_graph, t, s)
    
    parentsS = graphBFS["parent"]
    distS = graphBFS["dist"]
    parentsT = reverseGraphBFS["parent"]
    distT = reverseGraphBFS["dist"]
    shortestLength = graphBFS["shortestPathDistance"]
    print("shortest_path_length is ", shortestLength)

    # test all vertices except s and t to be in the subgraph
    vertices_to_test = set(graph.keys()) - set([s, t]) 
    print("dist S", distS)
    print("dist T",  distT)

    shortest_paths_subgraph = defaultdict(set)
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
                 
                add_path_to_graph(path_S_to_V_to_T, shortest_paths_subgraph)
                shortest_path_vertices.update(path_S_to_V_to_T)

                # TODO: you can optimize here by subtracting vertices you added from vertices to test. 
                             
    print("shortest path subgraph is: ", shortest_paths_subgraph)

    return {
        "shortest_path_dag": shortest_paths_subgraph,
        "shortest_path_vertices": shortest_path_vertices, 
        "shortest_path": shortestLength
    }



def crazyBFS(graph, Z, shortest_paths_dag, vertices_in_dag):
    seen, queue = set([s]), deque([s])
    parent, dist = {}, {}
    dist[s] = 0
    parent[s] = None

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


def create_longer_path_using_lost_edges(graph, shortest_path_dag, vertices_in_dag):
    for V in vertices_in_dag: 
        all_neighbors = graph[v]
        dag_neighbors = shortest_path_dag[v]

        not_in_dag_neighbors = all_neighbors - dag_neighbors # set subtraction
        
        for K in not_in_dag_neighbors:
            if(k in vertices_in_dag):
                # We found a lost edge. REPORT THAT 2 PATHS OF DIFFERENT SIZE EXIST.
                # Lost edge creates [S->V, V->K, K->T] which is longer than shortest path 
                # Question why would K -> T exist? it has to exist because K is a vertice in the DAG
                return True
    
    return False

    
def solution(graph, s, t):

    # Reversed graph is needed. 
    reversed_graph = defaultdict(set) # a map with values as list type

    # Go through graph and create reversed graph
    for parent, neighbors in graph.items():
        for n in neighbors:
            reversed_graph[n].add(parent)
    
    pprint.pprint(reversed_graph)

    buildResult = create_shortest_paths_subgraph(graph, reversed_graph, s, t)
    shortest_path_dag = buildResult["shortest_path_dag"] # subgraph that contains all shortest paths from s to t
    shortest_path_dag_vertices = buildResult["shortest_path_vertices"]
    shotest_path = buildResult["shortest_path"] # shortest path from s to t

    




# DOES THERE EXIST AT LEAST 2 DIRECTED SIMPLE PATHS FROM S TO T of different lengths
solution(graph=g, s=1, t=5)
    



