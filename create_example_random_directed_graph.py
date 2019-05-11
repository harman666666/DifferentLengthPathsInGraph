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
        
        # the graph cannot have a node that has itself as a neighbor
        possible_neighbors = set(range(vertices)) - set([i])
        
        ne = np.random.randint(vertices) # degree is random

        g[i] = set(np.random.choice(list(possible_neighbors), ne, replace=False))


    return g

