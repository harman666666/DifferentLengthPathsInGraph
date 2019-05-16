import pprint

from utility import verify_path_exists
from brute_force_dfs_solution import brute_force_solution
from polynomial_solution import poly_solution








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

print("################################# BRUTE FORCE SOLUTION")


brute_force_soln = brute_force_solution(g, 1, 5)

if(brute_force_soln["result"]):
    a = verify_path_exists(g, brute_force_soln["a_shortest_path"])
    b = verify_path_exists(g, brute_force_soln["a_longer_path"])
    if(a):
        print("brute force shorter path was correct")
    else: 
        print("brute force shorter path was NOTTTT correct")
    
    if(b):
        print("brute force longer path was correct")
    else: 
        print("brute force longer path was NOTTTT correct")


else:
    print("Solution was not found with brute force solution")

print("################################# POLY SOLUTION")

poly_soln = poly_solution(g, 1, 5)

if(poly_soln["result"]):
    c = verify_path_exists(g, poly_soln["a_shortest_path"])
    d = verify_path_exists(g, poly_soln["a_longer_path"])
    if(c):
        print("poly solution shorter path was correct")
    else: 
        print("poly solution shorter path was NOTTTT correct")
    
    if(d):
        print("poly solution longer path was correct")
    else: 
        print("poly solution longer path was NOTTTT correct")


else:
    print("Solution was not found with poly solution")