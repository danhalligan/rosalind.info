# Shortest Cycle Through a Given Edge

from .helpers import parse_graphs
from .dij import dij

# For this problem, we start Dijaktra's search at the end of the specified
# edge, and take the distance to the start of the edge (if reachable) and add
# the weight of the specified edge to get the cycle length.


def first_edges(handle):
    lines = handle.read().splitlines()
    ngraphs = int(lines[0])
    start = 1
    edges = []
    for i in range(ngraphs):
        if lines[start] == "":
            start += 1
        _, n_edges = lines[start].split()
        edges += [list(map(int, lines[start + 1].split()))]
        start += int(n_edges) + 1
    return edges


def cte(graph, edge):
    dist = dij(graph, start=edge[1])[edge[0] - 1]
    return dist if dist == -1 else dist + edge[2]


def main(file):
    edges = first_edges(open(file))
    graphs = list(parse_graphs(open(file), directed=True, weighted=True))
    res = [cte(graphs[i], edges[i]) for i in range(len(graphs))]
    print(*res)
