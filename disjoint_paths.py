


'''

Find path from s to t



'''

from collections import deque
from utility import bfs_with_restriction_set, \
                    dfs_with_restriction_set, \
                    get_path_to_root

def brute_force_disjoint_paths(graph, s, t, y, z, DEBUG=False):
    # CREATE path from s to t
    
    # other algo tries to create path from y to z
    # if it sees vertex from s to t, 
    # saves the location right, and tries to complete dfs,
    # if it fails ->

    '''
    go back to saved locations 
    WHERE WE INTERSECTING IN OUR ORIGINAL
    DFS,
    s -> t bfs, we tell it to choose a different neighbor, then BFS to t.
    then we try dfs from y to z again starting at that node we just changed!
    
    
    
    
    if it still didnt work, go to 
        
    '''


    s_to_t_bfs = bfs_with_restriction_set(graph, s, t, set())

    if(s_to_t_bfs["result"] == False):
        if(DEBUG): print("FAILED ON FIRST DFS. from s -> x", (s, X))
        return {
            "result": False
        }

    
    s_to_t_path = get_path_to_root(s_to_t_bfs["parents"], t)[::-1]
    '''
    # TRY TO DFS FROM Y TO Z!
    # We may touch a vertex in s->t path
    # if we fail, go back to these touch points,

    # recreate path from s to t, 
    # go through each touch point, and make path that avoids touch point!
    # from s to t, using bfs.

    # try dfs again, we will get new touch points!
    # add new touch points to old touch points, and modify s -> t path,
    # and try to make path from y to z again!
    
    dfs_result = dfs_with_restriction_set(graph, s, t, set(s_t_path))
    if dfs_result == True => RETURN SOLUTION
    else:
        # find places where dfs touchd!
        touches = stack( zip(dfs_result.touches(), [array full of st-bfs path]) )
        checked_touch_point = set()
        for (touch, old_st_path) in touches:
            
            if(touch in checked_touch_point):
                continue
            else:
                checked_touch_point.add(touch)
            
            split bfs s_t path into  2 pieces, part before touch, 
                                                and part after touch

            new_s_t_path = => node before touch point (call it Z)
                           => choose different neighbor and bfs to t
                           => if failed, move on to next neighbor
                           => in list of neighbors for Z
                           => all neighbors exhausted, 
                                move on to next touch point.
            
            NEW_TOUCH_POINTS = TRY TO DFS AGAIN FROM Y TO Z
            append new touch points to stack!
            stack.append(zip(touch_points, old_st_path + choosen_neighbor ))

            IF WORKS => RETURN TRUE
            ELSE IF FAIL, CONTINUE WITH NEXT TOUCH POINT!



    for touchpoint in touch_points_from
    
    '''


    # ALGORITHM 2:

    '''

    bfs from s to t

    bfs from y to z 

    BFS does shortest path, so each vertex in bfs is very important

    okay so common vertices between 2 bfs paths are called critical vertices!

    one path can have a critical vertex, but then other path has to avoid 
    
    go through each common vertex, use as part of restriction set. 
    .... THIS MIGHT BE GOOD.    

    '''



    
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
   
print ("The maximum possible flow is %d " % g.FordFulkerson(source, sink)) 
  
#This code is contributed by Neelam Yadav 