import pprint

from utility import verify_solution_if_paths_exist, \
                    create_example_rand_directed_graph, \
                    create_graph_without_longer_path,  \
                    create_graph_with_h_path

from brute_force_dfs_solution import brute_force_solution
from polynomial_solution import poly_solution

import pickle
import time
import datetime
import os

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


def test_when_shortest_path_is_length_1():
    g =  {0: set([19, 2, 67, 70, 6, 40, 42, 99, 21, 54, 57]), 1: set([0, 97, 98, 99, 22, 51, 71, 72, 16, 50, 67, 52, 86, 89, 90, 63, 60, 69, 31]), 2: set([98, 67, 7, 44, 46, 48, 17, 82, 83, 78, 25, 95, 94, 63]), 3: set([42, 35, 37, 71, 41, 74, 11, 44, 79, 51, 53, 73, 88, 25, 61, 30]), 4: set([6, 7, 76, 45, 46, 59, 29]), 5: set([66, 22, 38, 41, 74, 48, 86, 73, 92, 30]), 6: set([65, 98, 99, 68, 73, 11, 13, 78, 15, 80, 17, 52, 85, 34, 55, 88, 58, 91]), 7: set([96, 66, 27, 37, 74, 47, 19, 52, 21, 91, 94]), 8: set([65, 99, 36, 39, 72, 74, 51, 84, 85, 22, 90, 29, 62, 63]), 9: set([]), 10: set([68, 70, 7, 73, 75, 12, 77, 19, 85, 87, 61, 62]), 11: set([0, 68, 38, 71, 72, 41, 74, 18, 16, 81, 82, 57, 91, 93, 31]), 12: set([0, 65, 67, 39, 44, 79, 49, 87, 89, 27, 60, 94]), 13: set([68, 71, 41, 63, 80, 81, 50, 56, 26, 91, 93, 94, 95]), 14: set([32, 97, 99, 36, 5, 39, 74, 11, 76, 13, 80, 19, 23, 24, 27, 10]), 15: set([57, 13, 39]), 16: set([66, 34, 36, 7, 40, 2, 77, 79, 48, 19, 20, 24, 52, 71, 60, 30, 63]), 17: set([51, 61, 93]), 18: set([0, 98, 35, 68, 5, 73, 43, 44, 78, 15, 59, 19, 52, 46, 24, 89, 58, 27, 28]), 19: set([64, 1, 67, 37, 97, 10, 71, 45, 14, 13, 18, 3, 20, 24, 89, 26]), 20: set([64, 65, 49, 97, 41, 43, 48, 17, 18, 86, 89, 63]), 21: set([65, 4, 69, 72, 42, 62, 78, 46, 15, 81, 82, 85, 23, 56, 89, 58, 63, 30, 31]), 22: set([32, 54, 38, 7, 8, 74, 43, 12, 15, 17, 82, 51, 52, 86, 24, 20, 60]), 23: set([32, 1, 67, 68, 5, 39, 10, 75, 49, 35, 30, 57, 59, 60, 94, 95]), 24: set([36]), 25: set([34, 36, 16, 56, 83, 54, 88, 27, 31]), 26: set([97, 34, 59, 68, 69, 33, 42, 12, 13, 78, 77, 27, 79, 56, 15, 31]), 27: set([0, 4, 9, 74, 51, 53, 60, 31]), 28: set([89, 5]), 29: set([33, 91]), 30: set([0, 83, 6, 87]), 31: set([16, 95, 7, 76, 63]), 32: set([]), 33: set([1, 2, 99, 6, 40, 74, 75, 44, 15, 49, 50, 35, 21, 22, 87, 56, 67, 58, 95]), 34: set([0, 88, 27, 81]), 35: set([]), 36: set([1, 66, 3, 69, 40, 73, 79, 12, 15, 90, 27, 84, 56, 58, 91, 93, 62]), 37: set([2, 67, 68, 49, 8, 41, 91, 92, 81, 82, 46, 54, 55, 56, 25, 59, 60]), 38: set([4, 7, 40, 45, 78, 47, 17, 84, 85, 54, 56, 89, 29]), 39: set([]), 40: set([95]), 41: set([3, 75, 11, 34]), 42: set([64, 33, 2, 68, 0, 77, 88, 29, 95]), 43: set([0, 69, 38, 40, 88, 20, 85, 56, 31]), 44: set([89, 90, 17]), 45: set([0, 51, 12, 42, 55, 76, 61, 48, 17, 44, 83, 24, 9, 88, 89, 59, 92, 10]), 46: set([41]), 47: set([97, 67, 36, 5, 6, 73, 75, 77, 78, 17, 51, 84, 53, 88, 68, 26, 91, 61]), 48: set([16, 17]), 49: set([35, 15]), 50: set([17]), 51: set([64, 4, 70, 8, 9, 43, 15, 18, 85, 87, 56, 68, 26, 79, 60, 95]), 52: set([61, 7, 9, 74, 43, 50, 47, 80, 81, 82, 20, 88, 57, 26, 59, 58, 31]), 53: set([18, 4, 6, 12, 13, 47, 82, 19, 20, 21, 22, 55, 57, 59, 89, 30]), 54: set([16, 34, 69, 70, 1, 10, 7, 77, 80, 86, 87, 56, 89, 58, 28, 26, 94]), 55: set([82, 67, 21, 77]), 56: set([0, 37, 11, 12, 13, 50, 57, 28]), 57: set([65, 67, 1, 5, 51, 97, 9, 42, 76, 79, 16, 19, 52, 62, 89, 38, 92, 29, 94]), 58: set([2, 11, 81, 19, 52, 89, 90, 92]), 59: set([0, 1, 34, 35, 70, 40, 73, 11, 2, 46, 15, 3, 52, 85, 25, 58]), 60: set([0, 33, 98, 11, 49, 40, 41, 71, 34, 80, 24, 50, 20, 9, 56, 23, 87, 61, 94]), 61: set([65, 66, 97, 6, 39, 9, 10, 34, 48, 17, 83, 41, 87, 24, 25, 28, 62]), 62: set([34, 5, 6, 15, 82, 20, 21, 38, 92]), 63: set([41, 60, 22, 28]), 64: set([96, 99, 36, 69, 6, 39, 47, 75, 50, 46, 61, 48, 27, 82, 59, 92, 93]), 65: set([63]), 66: set([96, 83, 36, 69, 39, 73, 14, 13, 46, 56, 19, 84, 78, 54, 24, 25]), 67: set([2, 43, 72, 12, 75, 76, 81, 3, 53, 22, 55, 24, 89, 26, 60, 90]), 68: set([64, 66, 67, 4, 6, 76, 77, 45, 48, 21, 22, 24, 89, 15, 28]), 69: set([65, 66, 6, 40, 41, 10, 12, 13, 78, 80, 49, 82, 29, 31]), 70: set([35, 36, 42, 12, 47, 81, 83, 53]), 71: set([36, 37, 81, 50, 19, 20, 30, 63]), 72: set([0, 64, 66, 4, 6, 40, 95, 77, 17, 82, 21, 89, 59, 94, 53]), 73: set([98, 3, 68, 5, 72, 28, 42, 80, 53, 57, 88, 36, 27, 92, 85, 21]), 74: set([12]), 75: set([33, 98, 36, 11, 48, 51, 62, 89, 58, 27, 61, 30]), 76: set([96, 6, 77, 15, 80, 81, 21, 30]), 77: set([80, 43, 15, 48, 86, 55, 57, 26, 59]), 78: set([33, 66, 39, 11, 76, 14, 17, 19, 20, 54, 23, 27]), 79: set([35, 67, 99, 81, 39, 42, 76, 18, 93, 80, 49, 82, 51, 21, 55, 24, 90, 29]), 80: set([0, 65, 66, 67, 38, 49, 7, 56, 75, 22, 55, 88, 25]), 81: set([64, 1, 66, 75, 5, 38, 71, 8, 43, 17, 52, 53, 72]), 82: set([]), 83: set([98, 66, 5, 39, 8, 77, 88, 21, 24, 57, 91, 30]), 84: set([64, 48, 98, 3, 37, 7, 8, 73, 71, 45, 78, 11, 82, 21, 60, 69, 63]), 85: set([64, 16, 99, 6, 97, 47, 48, 18, 21, 23, 24, 58, 60, 90]), 86: set([38, 73, 11, 13, 48, 83, 21, 87, 26, 27, 63]), 87: set([27]), 88: set([36, 5, 46, 51, 4, 31, 62, 37]), 89: set([66, 36, 81, 8, 74, 75, 76, 77, 13, 72, 18, 86, 24, 68, 90]), 90: set([96, 24, 26, 92, 5]), 91: set([41, 11, 14, 95]), 92: set([0, 64, 76]), 93: set([]), 94: set([75, 44, 13, 78, 17, 50, 20, 22, 28, 62]), 95: set([33, 38, 11, 78, 84, 87, 20, 58, 28]), 96: set([19, 1, 98, 36, 68, 9, 42, 47, 50, 51, 54, 4, 58, 95]), 97: set([56]), 98: set([97, 34, 67, 4, 70, 74, 2, 46, 79, 49, 52, 53, 22, 57, 56, 88, 27, 94]), 99: set([4, 39, 8, 28, 10, 77, 81, 82, 85, 56, 90, 60, 26, 94, 63])}
    S = 1
    T = 89

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
    
    print("##############################################")
    print(brute_force_soln)
        

def create_performance_graphs(name, amount=1, vertices=10000, max_neighbors=100):
    graphs = []
    for i in range(amount):
        g = create_example_rand_directed_graph(vertices, max_neighbors)
        graphs.append(g)

    with open(name, 'wb') as f:
        pickle.dump(graphs, f)


def performance_testing(name, s=1, t=19):
    # One way is to create lots of big graphs, and test both algos and see which ones better
    with open(name, 'rb') as f:
        graphs = pickle.load(f)

    start_time = time.time()

    for g in graphs:
        run_poly_example(g, s, t)
    
    duration = (time.time() - start_time)
    print("TIME TO EXECUTE IN SECONDS FOR POLY WAS", duration)
    
    start_time2 = time.time()
    for g in graphs:
            run_brute_force_example(g, s, t)
    
    duration2 = (time.time() - start_time2)
    print("TIME TO EXECUTE IN SECONDS FOR BRUTE FORCE WAS", duration2)

def performance_testing_2():
    square_size = 100
    g = create_graph_without_longer_path(square_size)
    s = 0
    t = square_size** 2 - 1

    # run_example(g, S, T)

    start_time = time.time()

    run_poly_example(g, s, t)
    
    duration = (time.time() - start_time)
    print("TIME TO EXECUTE IN SECONDS FOR POLY WAS", duration)
    
    start_time2 = time.time()
    
    run_brute_force_example(g, s, t)
    
    duration2 = (time.time() - start_time2)
    print("TIME TO EXECUTE IN SECONDS FOR BRUTE FORCE WAS", duration2)





# SAVE RESULTS FOR THE TEST USING BRUTEFORCE!

def benchmark_correctness_testing(DEBUG=False):
    score = 0
    i = 0
    errors = 0
    # Create file that saves examples that lead to errors
    filename = "correctness-benchmarking/correctness-testing-" + \
                    str(datetime.datetime.now().date()) + '_' + \
                    str(datetime.datetime.now().time()).replace(':', '.')
    
    os.makedirs(os.path.dirname(filename), exist_ok=True)


    fp = open( filename, "w+")

    while True:
                
        #g = create_example_rand_directed_graph(vertices=50, max_neighbors=3)
        g = create_graph_with_h_path(20, 2)["g"]

        Sarr = [1, 2,3,4,5]
        Tarr = [10,11,12,13,14]

        # ITERATE THROUGH DIFFERENT S AND T:
        fail = False 
        for S in Sarr:
            for T in Tarr: 
                i += 1
                
                if DEBUG: pprint.pprint(g)

                if DEBUG: print("################################# BRUTE FORCE SOLUTION, index is ", i)


                brute_force_soln = brute_force_solution(g, S, T)
                if DEBUG:  print(brute_force_soln)
                if(brute_force_soln["result"]):
                    a = verify_solution_if_paths_exist(g, 
                                                       brute_force_soln["a_shortest_path"], 
                                                       brute_force_soln["a_longer_path"], S, T)


                    if(a):
                        if DEBUG: print("brute force correct")
                    else:
                        print("solution brute force came up with is WRONG")
                        break
                else:
                    if DEBUG: print("Solution was not found with brute force solution")

                if DEBUG: print("################################# POLY SOLUTION, index is ", i)

                poly_soln = poly_solution(g, S, T)

                if(poly_soln["result"]):
                    b = verify_solution_if_paths_exist(g, poly_soln["a_shortest_path"],poly_soln["a_longer_path"], S, T)
                    if(b):
                        if DEBUG: print("poly soln is correct")
                    else:
                        print("solution POLY came up with is WRONG")
                        break


                else:
                    if DEBUG: print("Solution was not found with poly solution")
                
                if DEBUG: print("##############################################")
                if DEBUG: print(brute_force_soln)
                
                if brute_force_soln["result"] != poly_soln["result"]:
                    print("BRUTE FORCE SOLUTION RESULT AND POLY SOLUTION RESULT DIFFER. BAD BREAK", file=fp)
                    print("BRUTE FORCE SOLUTION RESULT AND POLY SOLUTION RESULT DIFFER. BAD BREAK")

                    pprint.pprint(g, fp)
                    pprint.pprint(g)
                    print("S AND T WERE ", (S, T), file=fp)
                    print("S AND T WERE ", (S, T))

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

def run_poly_example(g, S, T):

    poly_soln = poly_solution(g, S, T)

    #if(poly_soln["result"]):
    #    b = verify_solution_if_paths_exist(g, poly_soln["a_shortest_path"],poly_soln["a_longer_path"], S, T)
    #    if(b):
    #        print("poly soln is correct")
    #    else:
    #        print("solution POLY came up with is WRONG")


    #else:
    #    print("Solution was not found with poly solution")

def run_brute_force_example(g, S, T):
    brute_force_soln = brute_force_solution(g, S, T)


def run_example(g, S, T):
    print("GRAPH IS: ")
    pprint.pprint(g)

    print("################################# BRUTE FORCE SOLUTION")


    brute_force_soln = brute_force_solution(g, S, T)

    if(brute_force_soln["result"]):
        a = verify_solution_if_paths_exist(g, brute_force_soln["a_shortest_path"], brute_force_soln["a_longer_path"], S, T)


        if(a):
            print("brute force correct")
            print(brute_force_soln)
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
  


# HARD EXAMPLE 1 ONLY PASSES IF YOU USE LOST EDGES METHOD
def hard_example_1():
    g = {0: set([]), 1: set([13, 14]), 2: set([14, 7]), 3: set([]), 4: set([]), 5: set([19, 4]), 6: set([11]), 7: set([10, 12]), 8: set([17, 11]), 9: set([10]), 10: set([17]), 11: set([]), 12: set([16, 19]), 13: set([4]), 14: set([10, 7]), 15: set([]), 16: set([]), 17: set([19]), 18: set([9, 4]), 19: set([16])}
    S = 1
    T = 19
    run_example(g, S, T)



def hard_example_2():
    g = {0: set([]), 1: set([0, 5]), 2: set([]), 3: set([19, 15]), 4: set([16]), 5: set([18, 10]), 6: set([]), 7: set([13, 14]), 8: set([7]), 9: set([17]), 10: set([9, 12]), 11: set([9, 12]), 12: set([17, 13]), 13: set([]), 14: set([]), 15: set([3, 6]), 16: set([]), 17: set([3, 15]), 18: set([16, 12]), 19: set([4])}
    S = 1
    T = 19
    run_example(g, S, T)

def hard_example_3():
    g = {0: set([16, 9]), 1: set([5]), 2: set([]), 3: set([16, 0, 21]), 4: set([]), 5: set([24, 1, 18, 8]), 6: set([10, 20]), 7: set([17, 2, 18, 13]), 8: set([24, 9, 10, 7]), 9: set([1, 11]), 10: set([17, 20]), 11: set([12]), 12: set([]), 13: set([24, 12, 5]), 14: set([20]), 15: set([9, 12]), 16: set([21]), 17: set([24, 19, 10, 3]), 18: set([24, 10, 7]), 19: set([16, 10, 4, 7]), 20: set([24]), 21: set([23]), 22: set([7, 15]), 23: set([]), 24: set([23])}
    S = 1
    T = 19
    run_example(g, S, T)

# Tests lost edges method
# This hard test will remind you that you can traverse not only one but multiple
# lost edges to create a longer path
# Allow multiple lost edges to be traversed in both the lost edges method and the outer vertex method. 

def hard_example_4():
    g = {0: set([12]), 1: set([20, 5, 13]), 2: set([21]), 3: set([19, 20]), 4: set([1, 10, 2, 17]), 5: set([0, 11, 4]), 6: set([]), 7: set([24, 11, 21]), 8: set([10, 12, 13]), 9: set([]), 10: set([8, 5, 23]), 11: set([17]), 12: set([]), 13: set([8, 9, 7]), 14: set([12, 5]), 15: set([16, 8, 2, 3]), 16: set([4, 12, 13, 14]), 17: set([]), 18: set([]), 19: set([3]), 20: set([12, 4]), 21: set([]), 22: set([]), 23: set([19, 11, 22, 7]), 24: set([12, 23])}
    S = 1
    T = 19

    run_example(g, S, T)

# Outer vertex test
def hard_example_5():
    g = {0: set([3]), 1: set([15, 4, 6, 7]), 2: set([0, 16, 14, 23]), 3: set([]), 4: set([]), 5: set([12]), 6: set([0, 10, 4, 18]), 7: set([9, 13]), 8: set([20]), 9: set([17, 2, 5]), 10: set([]), 11: set([24]), 12: set([20, 22]), 13: set([3]), 14: set([6, 23]), 15: set([0]), 16: set([6, 22, 7]), 17: set([16, 0, 21, 9]), 18: set([8, 15]), 19: set([24, 11, 21]), 20: set([10, 11, 15]), 21: set([17, 13]), 22: set([8, 1, 19, 16]), 23: set([2, 20, 22]), 24: set([4, 13])}
    S = 1 
    T = 19
    run_example(g, S, T)

def hard_example_6():
    g =  {0: set([17]), 1: set([17]), 2: set([9, 14]), 3: set([17, 11, 15]), 4: set([9, 7, 13, 17]), 5: set([12, 13]), 6: set([20, 21]), 7: set([8]), 8: set([9, 3, 21, 14]), 9: set([16, 17, 18, 4]), 10: set([8, 24, 3, 7]), 11: set([0, 9, 23]), 12: set([3, 5, 6]), 13: set([16, 9, 19, 15]), 14: set([8, 9, 18, 7]), 15: set([9, 3, 14]), 16: set([3]), 17: set([16]), 18: set([]), 19: set([]), 20: set([9]), 21: set([7]), 22: set([9, 14, 15]), 23: set([3, 4]), 24: set([16, 11, 7])}
    S = 1 
    T = 19
    run_example(g, S, T)


def hard_example_7():
    g = {0: set([7]), 1: set([9, 10]), 2: set([1]), 3: set([17]), 4: set([5]), 5: set([11]), 6: set([19, 14]), 7: set([11, 19]), 8: set([9]), 9: set([2, 15]), 10: set([16, 11]), 11: set([0, 12]), 12: set([4, 6]), 13: set([1, 6]), 14: set([3]), 15: set([12]), 16: set([13]), 17: set([18, 19]), 18: set([8, 1]), 19: set([2, 10])}
    S = 1
    T = 14
    run_example(g, S, T)


def hard_example_8():
    '''
    BRUTE FORCE SOLUTION RESULT AND POLY SOLUTION RESULT DIFFER. BAD BREAK
    defaultdict(<class 'set'>,
            {0: {1, 10},
             1: {9},
             2: {8, 3},
             3: {9, 11},
             4: {3},
             5: {0},
             6: {10, 11},
             7: {2},
             8: {5},
             9: {0, 4},
             10: {7},
             11: {6}})
    S AND T WERE  (2, 10)
    '''
    # CREATE THIS!
    pass

benchmark_correctness_testing()
 
# performance_testing_2()

#test_when_shortest_path_is_length_1()
#hard_example_1()
#hard_example_2()

#hard_example_3()

# create_performance_graphs("1graph10000V100N")
# performance_testing(name="1graph10000V100N")

#hard_example_4()

# Outer vertex test
# hard_example_5()
# hard_example_6()

# hard_example_7()
