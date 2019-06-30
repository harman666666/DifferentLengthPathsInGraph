
from brute_force_get_all_paths import brute_force_get_all_solution


'''

Find path from s to t



'''

from collections import deque
from utility import bfs_with_restriction_set, \
                    dfs_with_restriction_set, \
                    get_path_to_root, \
                    create_example_rand_directed_graph, \
                    create_graph_without_longer_path,  \
                    create_graph_with_h_path

def brute_force_disjoint_paths(graph, s, t, y, z, DEBUG=False):

    all_paths_from_s_to_t = brute_force_get_all_solution(graph, s, t, get_shortest_paths_too=True)
    all_paths_from_y_to_z = brute_force_get_all_solution(graph, y, z, get_shortest_paths_too=True)

    print("all_paths_from_s_to_t", all_paths_from_s_to_t)
    print("all paths from y to z", all_paths_from_y_to_z)

    for i in all_paths_from_s_to_t["longer_paths"]:
        for j in all_paths_from_y_to_z["longer_paths"]:
            # print("i", i)
            # print("j", j)
            
            path_i = set(i)
            path_j = set(j)

            intersection = path_i.intersection(path_j)
            # print("intersection", intersection)
            if(len(intersection) == 0):
                return {
                    "is_disjoint": True,
                    "paths": [i, j] 
                }
    return {
        "is_disjoint": False,
        "paths": []
    }
            






g = create_example_rand_directed_graph(vertices=20, max_neighbors=4)

result = brute_force_disjoint_paths(g, 0, 10, 1, 19)


print(result)

    
'''

FOR MORE EFFICIENT TRY MAX FLOW, OR TRY AVOIDANCE STRATEGY FOR SEARCHING!
LETS TRY MAX FLOW!!1







'''

def disjoint_graph(graph, start1, end1, start2, end2):
    pass







# Python program for implementation of Ford Fulkerson algorithm 
   
from collections import defaultdict 
   
#This class represents a directed graph using adjacency matrix representation 
class Graph: 
   
    def __init__(self,graph): 
        self.graph = graph # residual graph 
        self. ROW = len(graph) 
        #self.COL = len(gr[0]) 
          
   
    '''Returns true if there is a path from source 's' to sink 't' in 
    residual graph. Also fills parent[] to store the path '''
    def BFS(self,s, t, parent): 
  
        # Mark all the vertices as not visited 
        visited =[False]*(self.ROW) 
          
        # Create a queue for BFS 
        queue=[] 
          
        # Mark the source node as visited and enqueue it 
        queue.append(s) 
        visited[s] = True
           
         # Standard BFS Loop 
        while queue: 
  
            #Dequeue a vertex from queue and print it 
            u = queue.pop(0) 
          
            # Get all adjacent vertices of the dequeued vertex u 
            # If a adjacent has not been visited, then mark it 
            # visited and enqueue it 
            for ind, val in enumerate(self.graph[u]): 
                if visited[ind] == False and val > 0 : 
                    queue.append(ind) 
                    visited[ind] = True
                    parent[ind] = u 
  
        # If we reached sink in BFS starting from source, then return 
        # true, else false 
        return True if visited[t] else False
              
      
    # Returns tne maximum flow from s to t in the given graph 
    def FordFulkerson(self, source, sink): 
  
        # This array is filled by BFS and to store path 
        parent = [-1]*(self.ROW) 
  
        max_flow = 0 # There is no flow initially 
  
        # Augment the flow while there is path from source to sink 
        while self.BFS(source, sink, parent) : 
  
            # Find minimum residual capacity of the edges along the 
            # path filled by BFS. Or we can say find the maximum flow 
            # through the path found. 
            path_flow = float("Inf") 
            s = sink 
            while(s !=  source): 
                path_flow = min (path_flow, self.graph[parent[s]][s]) 
                s = parent[s] 
  
            # Add path flow to overall flow 
            max_flow +=  path_flow 
  
            # update residual capacities of the edges and reverse edges 
            # along the path 
            v = sink 
            while(v !=  source): 
                u = parent[v] 
                self.graph[u][v] -= path_flow 
                self.graph[v][u] += path_flow 
                v = parent[v] 
  
        return max_flow 
  
   
# Create a graph given in the above diagram 
  
graph = [[0, 16, 13, 0, 0, 0], 
        [0, 0, 10, 12, 0, 0], 
        [0, 4, 0, 0, 14, 0], 
        [0, 0, 9, 0, 0, 20], 
        [0, 0, 0, 7, 0, 4], 
        [0, 0, 0, 0, 0, 0]] 
  
g = Graph(graph) 
  
source = 0; sink = 5
   
# print ("The maximum possible flow is %d " % g.FordFulkerson(source, sink)) 
  
#This code is contributed by Neelam Yadav 