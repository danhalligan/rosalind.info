# Find an Eulerian Cycle in a Graph

from copy import deepcopy


def find_cycle(graph, key):
    cycle = []
    cycle += [key]
    while len(graph[key]):
        key = graph[key].pop(0)
        cycle += [key]
    for k in list(graph):
        if len(graph[k]) == 0:
            del graph[k]
    return cycle


def eulerian_cycle(graph):
    """
    Find an Eulerian cycle (if it exists).
    """
    g = deepcopy(graph)
    cycle = find_cycle(g, list(g.keys())[0])
    while len(g):
        key = [x for x in cycle if x in g][0]
        i = cycle.index(key)
        cycle = cycle[:i] + find_cycle(g, key) + cycle[(i + 1) :]
    return cycle


def parse_edges(edges):
    """
    Parse adjacency list to graph
    """
    g = {}
    for edge in edges:
        k, v = edge.split(" -> ")
        g[k] = v.split(",")
    return g


def main(file):
    edges = open(file).read().splitlines()
    cycle = eulerian_cycle(parse_edges(edges))
    print("->".join(cycle))
