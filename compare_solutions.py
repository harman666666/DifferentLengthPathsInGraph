import pprint

from utility import verify_solution_if_paths_exist
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
    a = verify_solution_if_paths_exist(g, brute_force_soln["a_shortest_path"], brute_force_soln["a_longer_path"], 1, 5)


    if(a):
        print("brute force correct")
    else:
        print("solution brute force came up with is WRONG")
else:
    print("Solution was not found with brute force solution")

print("################################# POLY SOLUTION")

poly_soln = poly_solution(g, 1, 5)

if(poly_soln["result"]):
    b = verify_solution_if_paths_exist(g, poly_soln["a_shortest_path"],poly_soln["a_longer_path"], 1, 5)
    if(b):
        print("poly soln is correct")
    else:
        print("solution POLY came up with is WRONG")


else:
    print("Solution was not found with poly solution")