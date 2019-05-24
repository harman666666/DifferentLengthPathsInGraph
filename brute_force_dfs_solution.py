
from utility import bfs, get_path_to_root



# Does not work for cycles. 
def find_longer_path(graph, s, t, shortest_length):

    if s == t:
        return {
            "found_longer_path": False, 
        }
    
    
    stack = [ (s, [s], set([s]) ) ]

    while stack:
        arg = stack.pop()

        node = arg[0]
        curr_path = arg[1]
        seen = arg[2]

        neighbors = graph[node]
        for neighbor in neighbors:
            if(neighbor == t):
                if len(curr_path) > shortest_length:
                # Found longer path!
                    return {
                    "found_longer_path": True,
                    "a_longer_path": curr_path + [neighbor]
                    }
                # ignore if we just found the shortest path
                continue
            
            if(neighbor not in seen):
                newSeen =  set( list(seen) + [neighbor]) # deep copy
                stack.append((neighbor, curr_path + [neighbor], newSeen))
            else: 
                # Detected cycle. kill recursion
                continue
    
    return {
        "found_longer_path": False
    }

def brute_force_solution(graph, s, t, DEBUG=False):

    graphBFS = bfs(graph, s)

    bfs_parents_root_s = graphBFS["parents"]
    bfs_dist_root_s = graphBFS["dist"]

    shortest_length = bfs_dist_root_s.get(t)

    if(shortest_length is None):
        # There is no path from s to t, FAIL!
        if(DEBUG): print("Hi brute force here. THERE IS NO SHORTEST PATH. FAIL.")
        return {
            "result": False,     
        }


    a_shortest_path = get_path_to_root(bfs_parents_root_s, t)[::-1]


    '''
    use dfs to find longer path!
    start at s, 
    '''

    result = find_longer_path(graph, s, t, shortest_length)
    if(result["found_longer_path"]):
        if(DEBUG): print("FOUND A SHORTER AND LONGER SIMPLE PATH")
        if(DEBUG):  print("THE SHORTEst PATH IS THE FOLLOWING: ", a_shortest_path)
        if(DEBUG): print("the longer path is the following: ", result["a_longer_path"])
        return {
            "result": True,
            "a_shortest_path": a_shortest_path, 
            "a_longer_path": result["a_longer_path"]
            }
    else:
        if(DEBUG): print("DID NOT FIND 2 SIMPLE PATHS")
        return {
            "result": False,
           
        }



