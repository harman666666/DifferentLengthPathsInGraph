

import numpy as np
import pprint
from collections import defaultdict, deque


# Creates random directed unweighted graph. 
def create_example_rand_directed_graph(vertices, max_neighbors):
    
    # This is how we do it with 
    # an adjacency matrix representation
    # adjacency = np.random.randint(0,2,(vertices,vertices))

    g = defaultdict(set) # default dict is a key-value map datatype with values that have set type.
    
    # Now create adjacency list representation from this. Map is adjacency list
    for i in range(vertices):
        
        # the graph cannot have a node that has itself as a neighbor
        possible_neighbors = set(range(vertices)) - set([i])
        
        ne = np.random.randint(max_neighbors) # degree is random

        g[i] = set(np.random.choice(list(possible_neighbors), ne, replace=False))


    return g


def bfs(graph, s):

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
        "parents": parent,
        "dist": dist
    })

def dfs_with_restriction_set(graph, start, end, restriction_set, seen = None, parents=None):
    # dont dfs in a cycle. visited set stops that!

    if(seen is None):
        seen = set([start])
        parents = {}
        parents[start] = None # start has no parents

    if(start in restriction_set): # or start in seen):
        # Kill this branch. 
        return {
            "result": False,
        }
    
    seen.add(start)

    
    if(start == end):
        return {
                "result": True,
                "parents": parents,
                "seen": seen, 
                "start": start, 
                "end" : end
            }

    neighbors = graph[start]
    
    for i in neighbors:
        if( i not in seen): 
            parents[i] = start

            dfs_result = dfs_with_restriction_set(graph, i, end, restriction_set, seen, parents)
            if(dfs_result["result"]):
                return dfs_result

    return {
        "result": False
    }

# Get the path to the root of a DFS or BFS tree. 
# bfs_tree is also the parents hashmap. we climb parents to get to root!
def get_path_to_root(bfs_tree, node, DEBUG=False):
    # root has parent called undefiend
    n = node 
    path = [n]

    if(DEBUG): print("GET PATH TO ROOT CALLED WITH FOLLOWING; ")
    if(DEBUG): print("NODE IS: ", node)
    if(DEBUG): print("bfs_tree", bfs_tree)
    

    while True:
        if(DEBUG): print("n is", n)
        parent = bfs_tree[n]
        #print(parent)    

        if(parent is None):
            break

        path.append(parent)
        n = parent


    return path

def verify_solution_if_paths_exist(graph, small_path, long_path, s, t):

    if(len(small_path) >= len(long_path)):
        print("One path is not longer than the other. FALSE SOLUTION")
        return False
    

    if(not verify_path_is_simple(small_path)):    
        print("Small path is not simple")
        return False
    
    if(not verify_path_is_simple(long_path)):
        print("Long path is not simple")
        return False

    a = verify_path_exists(graph, small_path, s, t)
    b = verify_path_exists(graph, long_path, s, t)
    
    

    if(a and b):
        print("Two paths are GOOD. VERIFIED")
        return True
    elif(a):
        print("only short path was correct")
        return False
    elif(b): 
        print("only long path was correct")
        return False
    else: 
        print("neither path was correct")
        return False

def verify_path_is_simple(path):
    return len(path) == len(set(path))

def verify_path_exists(graph, path, s, t):
    
    parent = path[0]

    for i in range(1, len(path)):
        child = path[i]

        children = graph[parent]
        if(child not in children):
            print("FAILED")
            return False
        
        parent = child

    if(path[0] == s and path[-1] == t):
        print("THE PATH " + str(path) + " has been verified")
        return True
    else: 
        print("The path does not start at S, and end at T")
        return False



