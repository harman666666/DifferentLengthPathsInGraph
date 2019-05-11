

import numpy as np
import pprint
from collections import defaultdict, deque


# Creates random directed unweighted graph. 
def create_example_rand_directed_graph(vertices):
 
    # This is how we do it with 
    # an adjacency matrix representation
    # adjacency = np.random.randint(0,2,(vertices,vertices))

    g = defaultdict(set) # default dict is a key-value map datatype with values that have set type.
    
    # Now create adjacency list representation from this. Map is adjacency list
    for i in range(vertices):
        
        # the graph cannot have a node that has itself as a neighbor
        possible_neighbors = set(range(vertices)) - set([i])
        
        ne = np.random.randint(vertices) # degree is random

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
        "parent": parent,
        "dist": dist
    })

# Get the path to the root of a DFS or BFS tree. 
# bfs_tree is also the parents hashmap. we climb parents to get to root!
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