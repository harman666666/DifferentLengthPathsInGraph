import pprint

from utility import verify_solution_if_paths_exist, create_example_rand_directed_graph
from brute_force_dfs_solution import brute_force_solution
from polynomial_solution import poly_solution


def test_outer_vertex_method():
    # This graph will test outer vertex method for the polynomial solution
    g = {
        0 : set([1,2]),
        1 : set([2, 3]),
        2: set([3,]),
        3: set([0, 5, 7]),
        5: set([9]),
        9 : set([]),
        7 : set([])
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
    
def test_lost_edge():
    g = {
    "R": set(["A", "Q"]),
    "Q": set(["A", "B"]),
    "A":  set(["B"]),
    "B": set(["C", "D", "E"]),
    "H": set(["Q", "A", "X"]),
    "E": set(["H"]),
    "G": set(["X"]),
    "F": set(["X"]),
    "D": set(["E", "G"]),
    "C": set(["F"]), 
    "X": set(["R"])   
    }

    S = "A"
    T = "X"

    print("GRAPH IS: ")
    pprint.pprint(g)

    print("################################# BRUTE FORCE SOLUTION")


    brute_force_soln = brute_force_solution(g, S, T)

    if(brute_force_soln["result"]):
        a = verify_solution_if_paths_exist(g, brute_force_soln["a_shortest_path"], brute_force_soln["a_longer_path"], S, T)


        if(a):
            print("brute force correct")
        else:
            print("solution brute force came up with is WRONG")
    else:
        print("Solution was not found with brute force solution")

    print("################################# POLY SOLUTION")

    poly_soln = poly_solution(g, S, T)

    if(poly_soln["result"]):
        b = verify_solution_if_paths_exist(g, poly_soln["a_shortest_path"],poly_soln["a_longer_path"], S, T)
        if(b):
            print("poly soln is correct")
        else:
            print("solution POLY came up with is WRONG")


    else:
        print("Solution was not found with poly solution")

# test_lost_edge()


# SAVE RESULTS FOR THE TEST USING BRUTEFORCE!

def benchmark_correctness_testing():
    for i in range(1):
        g = create_example_rand_directed_graph(vertices=100, max_neighbors=20)
        S = 1
        T = 89

        print("GRAPH IS: ")
        pprint.pprint(g)

        print("################################# BRUTE FORCE SOLUTION")


        brute_force_soln = brute_force_solution(g, S, T)
        print(brute_force_soln)
        if(brute_force_soln["result"]):
            a = verify_solution_if_paths_exist(g, brute_force_soln["a_shortest_path"], brute_force_soln["a_longer_path"], S, T)


            if(a):
                print("brute force correct")
            else:
                print("solution brute force came up with is WRONG")
        else:
            print("Solution was not found with brute force solution")

        print("################################# POLY SOLUTION")

        poly_soln = poly_solution(g, S, T)

        if(poly_soln["result"]):
            b = verify_solution_if_paths_exist(g, poly_soln["a_shortest_path"],poly_soln["a_longer_path"], S, T)
            if(b):
                print("poly soln is correct")
            else:
                print("solution POLY came up with is WRONG")


        else:
            print("Solution was not found with poly solution")
        
        print("##############################################")
        print(brute_force_soln)
        

benchmark_correctness_testing()


