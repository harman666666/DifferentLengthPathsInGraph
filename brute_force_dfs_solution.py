
import pprint
from utility import bfs, get_path_to_root



# Does not work for cycles. 
def find_longer_path(graph, s, t, shortest_length):

    if s == t:
        return {
            "founder_longer_path": False, 
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
            
            






def brute_force_solution(graph, s, t):

    graphBFS = bfs(graph, s)

    bfs_parents_root_s = graphBFS["parent"]
    bfs_dist_root_s = graphBFS["dist"]

    shortest_length = bfs_dist_root_s.get(t)

    if(shortest_length is None):
        # There is no path from s to t, FAIL!
        print("Hi brute force here. THERE IS NO SHORTEST PATH. FAIL.")
        return False

    a_shortest_path = get_path_to_root(bfs_parents_root_s, t)[::-1]


    '''
    use dfs to find longer path!
    start at s, 
    '''

    result = find_longer_path(graph, s, t, shortest_length)
    if(result["found_longer_path"]):
        print("FOUND A SHORTER AND LONGER SIMPLE PATH")
        print("THE SHORTEst PATH IS THE FOLLOWING: ", a_shortest_path)
        print("the longer path is the following: ", result["a_longer_path"])
        return True


g = {

    0 : set([1,2]),
    1 : set([2, 3]),
    2: set([3,]),
    3: set([0, 5]),
    5: set([9]),
    9 : set([])
}

print("GRAPH IS: ")
pprint.pprint(g)
brute_force_solution(g, 1, 5)