
from brute_force_get_all_paths import brute_force_get_all_solution
from collections import deque
from utility import bfs_with_restriction_set, \
                    dfs_with_restriction_set, \
                    get_path_to_root, \
                    create_example_rand_directed_graph, \
                    create_graph_without_longer_path,  \
                    create_graph_with_h_path
import pickle
import time
import datetime
import os
import pprint


def brute_force_disjoint_paths( s, t, y, z, graph, DEBUG=False):

    all_paths_from_s_to_t = brute_force_get_all_solution(graph, s, t, get_shortest_paths_too=True)
    all_paths_from_y_to_z = brute_force_get_all_solution(graph, y, z, get_shortest_paths_too=True)

    #print("all_paths_from_s_to_t", all_paths_from_s_to_t)
    #print("all paths from y to z", all_paths_from_y_to_z)

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
                    "result": True,
                    "paths": [i, j] 
                }
    return {
        "result": False,
        "paths": []
    }
            
    
'''

FOR MORE EFFICIENT TRY MAX FLOW, OR TRY AVOIDANCE STRATEGY FOR SEARCHING!
LETS TRY MAX FLOW!!1


'''


def simple_disjoint_graph(s, t, y, z, graph, DEBUG=False):


    S_to_T_bfs = bfs_with_restriction_set(graph=graph, 
                                          start=s, 
                                          end=t, 
                                          restriction_set=set([y,z]))

    if(S_to_T_bfs["result"] == False):
        if(DEBUG): print("FAILED ON FIRST BFS. from s -> t", (s, t))
        return {
            "result": False
        }
    
    S_to_T_bfs_path = get_path_to_root(S_to_T_bfs["parents"], t)[::-1]
 
    # Do second dfs
    if DEBUG: print("S_to_T_bfs_path_found: ", S_to_T_bfs_path)
    
    Y_to_Z_dfs = dfs_with_restriction_set(graph=graph, 
                                          start=y, 
                                          end=z, 
                                          # avoidance_set=set(S_to_T_bfs_path),
                                          restriction_set=set(S_to_T_bfs_path))
    print("Y TO Z DFS RESULT IS:")
    print(Y_to_Z_dfs)

    if(Y_to_Z_dfs["result"] == True):
        if(DEBUG): print("FAILED ON FIRST DFS FROM Y TO Z")
        Y_to_Z_dfs_path = get_path_to_root(Y_to_Z_dfs["parents"], z)[::-1] 

        return {
               "result": True,
            "S_to_T_path": S_to_T_bfs_path,
            "Y_to_Z_path": Y_to_Z_dfs_path
        }

    Y_to_Z_bfs = bfs_with_restriction_set(graph=graph, 
                                          start=y, 
                                          end=z, 
                                          restriction_set=set([s,t]))

    if(Y_to_Z_bfs["result"] == False):
        if(DEBUG): print("FAILED ON SECOND BFS. from Y -> Z", (y, z))
        return {
            "result": False
        }

    Y_to_Z_bfs_path = get_path_to_root(Y_to_Z_bfs["parents"], z)[::-1]

    print("Y TO Z BFS", Y_to_Z_bfs_path)

    
    S_to_T_dfs_2 = dfs_with_restriction_set(graph=graph, 
                                            start=s, 
                                            end=t, 
                                            # avoidance_set=set(Y_to_Z_bfs_path),
                                            restriction_set=set(Y_to_Z_bfs_path))

    if(S_to_T_dfs_2["result"] == True):
        S_to_T_dfs_path = get_path_to_root(S_to_T_dfs_2["parents"], t)[::-1] 
        return {
               "result": True,
            "S_to_T_path": S_to_T_dfs_path,
            "Y_to_Z_path": Y_to_Z_bfs_path        
            }

    

    return {
        "result": False
    }



def disjoint_graph(s, t, y, z, graph, DEBUG=False):


    S_to_T_bfs = bfs_with_restriction_set(graph=graph, 
                                          start=s, 
                                          end=t, 
                                          restriction_set=set([y,z]))

    if(S_to_T_bfs["result"] == False):
        if(DEBUG): print("FAILED ON FIRST BFS. from s -> t", (s, t))
        return {
            "result": False
        }
    
    S_to_T_bfs_path = get_path_to_root(S_to_T_bfs["parents"], t)[::-1]
 
    # Do second dfs
    if DEBUG: print("S_to_T_bfs_path_found: ", S_to_T_bfs_path)
    
    Y_to_Z_dfs = dfs_with_restriction_set(graph=graph, 
                                          start=y, 
                                          end=z, 
                                          avoidance_set=set(S_to_T_bfs_path),
                                          restriction_set=set([s, t]))

    if(Y_to_Z_dfs["result"] == False):
        if(DEBUG): print("FAILED ON FIRST DFS FROM Y TO Z")
        return {
            "result": False,
        }

    Y_to_Z_dfs_path = get_path_to_root(Y_to_Z_dfs["parents"], z)[::-1] 
    
    if DEBUG: print("Y to Z dfs path found: ", Y_to_Z_dfs_path)

    # Do third dfs
    S_to_T_dfs = dfs_with_restriction_set(graph=graph, 
                                          start=s, 
                                          end=t, 
                                          avoidance_set=set(),
                                          restriction_set=set(Y_to_Z_dfs_path))
    
    if(S_to_T_dfs["result"] == True):
        if(DEBUG): print("IT WAS MERGED WITHIN 2 DFS'S")
        return {
            "result": True,
            "S_to_T_path": get_path_to_root(S_to_T_dfs["parents"], t)[::-1],
            "Y_to_Z_path": Y_to_Z_dfs_path
        }

    # The second dfs failed, so move on to the third dfs:
    if DEBUG: print("SECOND DFS FAILED. WE HAVE TO REPEAT BFS AGAIN. second was S->T", (s, t))
    ############################################################################
    #### TRY THE OTHER WAY -> 1 BFS, 2 DFS
   
    Y_to_Z_bfs = bfs_with_restriction_set(graph=graph, 
                                          start=y, 
                                          end=z, 
                                          restriction_set=set([s,t]))

    if(Y_to_Z_bfs["result"] == False):
        if(DEBUG): print("FAILED ON SECOND BFS. from Y -> Z", (y, z))
        return {
            "result": False
        }

    Y_to_Z_bfs_path = get_path_to_root(Y_to_Z_bfs["parents"], z)[::-1]
    
    if DEBUG: print("Y_to_Z_bfs_path_found: ", Y_to_Z_bfs_path)
    
    S_to_T_dfs_2 = dfs_with_restriction_set(graph=graph, 
                                            start=s, 
                                            end=t, 
                                            avoidance_set=set(Y_to_Z_bfs_path),
                                            restriction_set=set([y, z]))

    if(S_to_T_dfs_2["result"] == False):
        if(DEBUG): print("FAILED ON FIRST DFS FROM Y TO Z")
        return {
            "result": False,
        }

    S_to_T_dfs_2_path = get_path_to_root(S_to_T_dfs_2["parents"], t)[::-1] 
    
    if DEBUG: print("S_to_T_dfs_2_path: ", S_to_T_dfs_2_path)

    Y_to_Z_dfs_2 = dfs_with_restriction_set(graph=graph, 
                                            start=y, 
                                            end=z, 
                                            avoidance_set=set(),
                                            restriction_set=set(S_to_T_dfs_2_path))
    
    if(Y_to_Z_dfs_2["result"] == True):
        if(DEBUG): print("IT WAS MERGED WITHIN 4 DFS'S")
        return {
            "result": True,
            "S_to_T_path": S_to_T_dfs_2_path,
            "Y_to_Z_path": get_path_to_root(Y_to_Z_dfs_2["parents"], z)[::-1]
        }

    if DEBUG: print("FOURTH DFS FAILED. COMPLETE FAILURE", (y, z))

    return {
        "result": False
    }

def max_flow_disjoint_paths(s, t, y, z, graph, DEBUG=True):
    # create a residual graph!
    # each edge has weight 1!
    import copy

    residual = {}
    
    # We will be adding reverse edges to the residual graph as we do the max flow algorithm!


    '''
    Try this to get simple paths if ford fulk dont make em:

    Split each node v in the graph into to nodes: vin and vout.
    For each node v, add an edge of capacity one from vin to vout.
    Replace each other edge (u, v) in the graph with an edge from uout to vin of capacity 1.
    '''

    for i, neighbors in graph.items():

        i_out = str(i) + "o" # split vertex into 2. out and in version. 
        i_in = str(i) + "i" 
        residual[i_in] = { i_out: 1}
        residual[i_out] = {}

        for n in neighbors:
            n_in = str(n) + "i"
            if(residual.get(i_out) is None):
                residual[i_out] = { n_in : 1 }
            else: 
                residual[i_out][n_in] = 1 # VALUE IS CURRENT CAPACITY IN RESIDUAL GRAPH


    s = str(s) + "i"
    y = str(y) + "i"
    t = str(t) + "o"
    z = str(z) + "o"
    
    # Connect A to all start nodes

    # Connect B to all end nodes
    # Max flow from A to B, get all edge independent paths from A to B
    # then get min cut edges! == max flow

    # then check if you can use min cut edges to 
    # create independent paths from s to t, and y to z?
    residual["A"] = {}
    residual["B"] = {}
    residual["A"][s] = 1
    residual["A"][y] = 1
    residual[t]["B"] = 1
    residual[z]["B"] = 1 

    original_graph = copy.deepcopy(residual)


    # remove in and out str parts from path!
    def normalize_path(path):
        # removes the in and out part from the path!
        p = []
        for i in path:
            if i[-1] == "i":
                p.append(i[:-1])
            else:
                continue
        
        return p


    if DEBUG: pprint.pprint(residual)

    '''Returns true if there is a path from source 's' to sink 't' in 
    residual graph. Also fills parent[] to store the path '''
    def BFS(s, t): 
        from collections import deque

        # Mark all the vertices as not visited 
        visited = set()
        parent = {}
        # Create a queue for BFS 
        visited, queue = set([s]), deque([s]) 
        parent[s] = None

         # Standard BFS Loop 
        while queue: 
  
            #Dequeue a vertex from queue and print it 
            u = queue.pop() 
          
            # Get all adjacent vertices of the dequeued vertex u 
            # If a adjacent has not been visited, then mark it 
            # visited and enqueue it 
            for neighbor in residual[u]: 
                if neighbor not in visited and residual[u][neighbor] > 0 : 
                    queue.appendleft(neighbor) 
                    visited.add(neighbor) 
                    parent[neighbor] = u 
  
        # If we reached sink in BFS starting from source, then return 
        # true, else false 
        return (True, parent) if t in visited else (False, parent)

    # go from s to t, then go from y to z
    # keep alternating, until doing both yield False?
    max_flow = 0

    # max flow variables for multiple source, sink
    finding_st_right_now = True
    source = None
    sink = None
    
    st_path = None
    yz_path = None

    while True: 
        
        if(finding_st_right_now):
            # source = s 
            # sink = t 
            source = "A"
            sink = "B"
            finding_st_right_now = False 
        else:
            # source = y 
            # sink = z 
            source = "A"
            sink = "B"
            finding_st_right_now = True


        print("source and sink are", (source, sink))
        (result, bfs_tree) = BFS(source, sink)
      
        if(result == False):
            if DEBUG: print("COULD NOT FIND PATH. BREAK")

            return {
            "result": False
            }
            break

        # Add path flow to overall flow 
        max_flow += 1

        path = get_path_to_root(bfs_tree, sink)[::-1]
        
        if(source == s):
            st_path = path
        else:
            yz_path = path

        # update graph!
        
        # update residual capacities of the edges and reverse edges 
        # along the path 
        v = sink
        while(v != source):
            u = bfs_tree[v]
            residual[u][v] -= 1
            

            # Do not create a backward edge from an Input to an outputself
            # Input edges should not have backward edges because they can ONLY Lead out of one output. 
            
            # if backward edge isnt there, add it!
            if DEBUG: print("v[-1]", v[-1])
            

            # only outs can connect to ins on the reverse side. 
            reverse_residual_v = v
            reverse_residual_u = u 
            if(v[-1] == "i"):
                # connect v-out to u-in 
                # even if we are dealing with v-in
                reverse_residual_v = v[:-1] + "o"
                reverse_residual_u = u[:-1] + "i"

             
            if residual[reverse_residual_v].get(reverse_residual_u) is None:
                residual[reverse_residual_v][reverse_residual_u] = 0
                
            residual[reverse_residual_v][reverse_residual_u] += 1
            v = bfs_tree[v]
        
        if DEBUG: print("FLOW FOUND IS", path)
        if DEBUG: print("NORMALIZED FLOW", normalize_path(path))
        if DEBUG: print("RESIDUAL AFTER FINDING A PATH from", (source, sink))
        if DEBUG: pprint.pprint(residual)

        if(max_flow == 2):
            if DEBUG: print("FOUND BOTH PATHS! WE GUCCI")

            break

    
    return_output = {
            "result": True,
            # "S_to_T_path": normalize_path(st_path),
            # "Y_to_Z_path": normalize_path(yz_path) # get_path_to_root(Y_to_Z_dfs_2["parents"], z)[::-1]
    }

     # original graph has all edge relationships with postiive val. 
     # if they now 0, its a min cut edge
    for i, neighbors in original_graph.items():
        for j in neighbors:
            if(residual[i][j] == 0):
                print("MINIUM CUT EDGE IS ", (i, j))                 

    if DEBUG: pprint.pprint(return_output)

    return return_output



    # build graph. has to be simple paths, so convert
    # each vertex to an edge of length 1






    # Do max flow algo (iterate back and forth between 2 nodes!)
    # start from first node -> go to target 1. 
    # start from second node -> go to target 2
    # if you fail! stop!



    

    








def benchmark_correctness_testing(DEBUG=False):
    score = 0
    i = 0
    errors = 0
    # Create file that saves examples that lead to errors
    filename = "disjoint-testing/d-testing-" + \
                    str(datetime.datetime.now().date()) + '_' + \
                    str(datetime.datetime.now().time()).replace(':', '.')
    
    os.makedirs(os.path.dirname(filename), exist_ok=True)


    fp = open( filename, "w+")

    while True:
                
        #g = create_example_rand_directed_graph(vertices=50, max_neighbors=3)
        g = create_example_rand_directed_graph(vertices=8, max_neighbors=2)

        # Sarr = [1, 2,3,4,5]
        # Tarr = [11,12,13,14, 15]
        # Yarr = [6,7,8,9,10]
        # Zarr = [16,17,18,19,20]

        Sarr = [1, 2]
        Tarr = [5,6]
        Yarr = [3,4]
        Zarr = [7,8]

        # ITERATE THROUGH DIFFERENT S AND T:
        fail = False 
        for (S,Y) in zip(Sarr, Yarr):
            for (T,Z) in zip(Tarr, Zarr): 
                i += 1
                
                if DEBUG: pprint.pprint(g)

                if DEBUG: print("################################# BRUTE FORCE SOLUTION, index is ", i)


                brute_force_soln = brute_force_disjoint_paths(S, T, Y, Z, g, DEBUG)
                # poly_soln = simple_disjoint_graph(S, T, Y, Z, g, DEBUG)
                poly_soln = max_flow_disjoint_paths(S, T, Y, Z, g, DEBUG)
                if brute_force_soln["result"] != poly_soln["result"]:
                    print("BRUTE FORCE SOLUTION RESULT AND POLY SOLUTION RESULT DIFFER. BAD BREAK", file=fp)
                    print("BRUTE FORCE SOLUTION RESULT AND POLY SOLUTION RESULT DIFFER. BAD BREAK")

                    pprint.pprint(g, fp)
                    pprint.pprint(g)
                    print("S AND T AND Y AND Z WERE ", (S, T, Y, Z), file=fp)
                    print("S AND T Y AND Z  WERE ", (S, T, Y, Z))

                    # fail = True
                    errors += 1
                    break
                else:
                    score += 1
                
                print("index: " +  str(i) + " is correct. poly soln matches brute force." + "tot errors: " + str(errors))

            #if(fail): 
            #    break
        if(fail): 
                break
    print("THE SCORE OUR POLY SOLUTION RECIEVED IS ", score)



def run_example(g, S, T, Y, Z, DEBUG):
    brute_force_soln = brute_force_disjoint_paths(S, T, Y, Z, g, DEBUG)
    # poly_soln = simple_disjoint_graph(S, T, Y, Z, g, DEBUG)
    poly_soln = max_flow_disjoint_paths(S, T, Y, Z, g, DEBUG)

    if brute_force_soln["result"] != poly_soln["result"]:

        
        print("BRUTE FORCE SOLUTION RESULT AND POLY SOLUTION RESULT DIFFER. BAD BREAK")
        
        pprint.pprint(g)
        print("S AND T Y AND Z  WERE ", (S, T, Y, Z))
        # fail = True
        
        print("BRUTE FORCE RESULT")
        pprint.pprint(brute_force_soln)
        print("poly soln")
        pprint.pprint(poly_soln)


    else:   
        print("SAME SOLUTION")
        print("BRUTE FORCE RESULT")
        pprint.pprint(brute_force_soln)
        print("poly soln")
        pprint.pprint(poly_soln)


def ex_1():
    g = {0: {8},
             1: {0, 6},
             2: {1, 18},
             3: {0, 11, 13},
             4: {17, 19},
             5: {19, 4},
             6: {1},
             7: {16, 13},
             8: {12, 5, 6},
             9: {8},
             10: {9, 6, 15},
             11: {0, 9, 3},
             12: {3, 14},
             13: {19, 15},
             14: {10},
             15: {14, 7},
             16: {8, 1, 14},
             17: {10, 20},
             18: set(),
             19: {4, 13},
             20: {9, 6, 1}}
    s = 2
    t = 15
    y= 7
    z = 20
    run_example(g, s, t, y, z, True)

# ex_1()

def simple_ex_1():
    g = {0: {2, 19},
             1: {9, 10, 17},
             2: {16, 12, 4},
             3: set(),
             4: {1},
             5: {8},
             6: {16, 19},
             7: {19, 20, 13},
             8: {16, 17, 12},
             9: {3, 14, 7},
             10: {8, 11, 6},
             11: {8, 18},
             12: set(),
             13: {2, 3},
             14: {6},
             15: {11, 3, 6},
             16: {8},
             17: {20, 4, 5},
             18: {20},
             19: {1, 10, 9},
             20: {2, 4}}
    s = 4
    t = 12
    y = 9
    z = 17

    run_example(g, s, t, y, z, True)




def simple_ex_2():
    g =  {0: {18, 5},
             1: {15},
             2: set(),
             3: {12},
             4: {15},
             5: {12, 14},
             6: {8, 2, 13},
             7: {15},
             8: {10},
             9: {12},
             10: {8, 5, 6},
             11: set(),
             12: {5},
             13: {17, 18, 14},
             14: {13},
             15: {12, 5},
             16: {19, 15},
             17: {19, 12},
             18: set(),
             19: {16, 18},
             20: {10, 2, 19}}
    s = 4
    t = 12
    y = 9
    z = 17

    run_example(g, s, t, y, z, True)

def simple_ex_3():
    '''
    BRUTE FORCE SOLUTION RESULT AND POLY SOLUTION RESULT DIFFER. BAD BREAK
defaultdict(<class 'set'>,
            {0: {5},
             1: {4},
             2: {0},
             3: {5},
             4: {7},
             5: set(),
             6: set(),
             7: {3}})
S AND T AND Y AND Z WERE  (1, 5, 3, 7)
    '''

    g = {0: {5},
             1: {4},
             2: {0},
             3: {5},
             4: {7},
             5: set(),
             6: set(),
             7: {3}}

    s = 1
    t = 5
    y = 3
    z = 7
    run_example(g, s, t, y, z, True)

# simple_ex_2()


simple_ex_1()
# benchmark_correctness_testing()
# simple_ex_3()