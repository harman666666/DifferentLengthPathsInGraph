import pprint

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


brute_force_solution(g, 1, 5)



print("################################# POLY SOLUTION")

poly_solution(g, 1, 5)
