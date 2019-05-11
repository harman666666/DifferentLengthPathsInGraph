
from utility import bfs, get_path_to_root


def brute_force_solution(graph, s, t):

    bfs(graph, s)

    '''
    find shortest path with bfs. 

    use dfs to find longer path!
    start at s, 
    '''