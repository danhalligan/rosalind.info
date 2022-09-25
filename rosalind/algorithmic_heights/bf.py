# Bellman-Ford Algorithm

from .helpers import parse_graph
from math import inf, isinf


def nodes(graph):
    s = list(graph.keys())
    e = [y["n"] for v in graph.values() for y in v]
    return set(s) | set(e)


def bf(graph, start=1):
    d = {n: inf for n in nodes(graph)}
    d[start] = 0
    n = nodes(graph)
    for _ in range(len(n) - 1):
        for u, x in graph.items():
            for v in x:
                if d[u] + v["w"] < d[v["n"]]:
                    d[v["n"]] = d[u] + v["w"]
    return ["x" if isinf(v) else v for k, v in d.items()]


def main(file):
    print(*bf(parse_graph(open(file), directed=True, weighted=True)))
