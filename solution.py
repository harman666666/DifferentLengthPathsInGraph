
import numpy as np
import pprint
from collections import deque  

# Creates random directed unweighted graph. 
def create_rand_directed_graph(vertices):
 
    # This is how we do it with 
    # an adjacency matrix representation
    # adjacency = np.random.randint(0,2,(vertices,vertices))

    g = {}
    
    # Now create adjacency list representation from this. Map is adjacency list
    for i in range(vertices):
        ne = np.random.randint(vertices) 
        g[i] = np.random.choice(vertices, ne, replace=False)


    return g
g = create_rand_directed_graph(10)

pprint.pprint(g)


def bfs(graph, s, t):

    seen, queue = set([s]), deque([s])
    parent, dist = {}, {}
    dist[s] = 0
    
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




# Shortest paths subgraph contains all the shortest paths from node s to t.
# it is a DAG from s to t. 
def create_shortest_paths_subgraph(graph, s, t):
    # You can do this by doing BFS on every vertex. Then
    pass

    


def solution(graph):
    pass
    



