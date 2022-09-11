# Topological Sorting

from .helpers import parse_graph
from collections import defaultdict
from .dag import graph_nodes


def topological_sort(g):
    def find_sort(v, visited, stack):
        visited[v] = True
        for i in g[v]:
            if not visited[i]:
                find_sort(i, visited, stack)
        stack.append(v)

    visited = defaultdict(bool)
    stack = []
    for i in graph_nodes(g):
        if not visited[i]:
            find_sort(i, visited, stack)

    return stack[::-1]


def main(file):
    g = parse_graph(open(file), directed=True)
    print(*topological_sort(g))
