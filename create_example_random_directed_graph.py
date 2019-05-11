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
