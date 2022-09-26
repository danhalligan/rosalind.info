# Negative Weight Cycle

from .bf import nodes, nedges
from .helpers import parse_graphs


def nwc(graph):
    dist = {n: 10**20 for n in nodes(graph)}
    dist[1] = 0
    for _ in range(nedges(graph) - 1):
        for u, x in graph.items():
            for v in x:
                dist[v["n"]] = min(dist[u] + v["w"], dist[v["n"]])
    for u, x in graph.items():
        for v in x:
            if dist[u] + v["w"] < dist[v["n"]]:
                return 1
    return -1


def main(file):
    graphs = parse_graphs(open(file), directed=True, weighted=True)
    print(*[nwc(g) for g in graphs])
