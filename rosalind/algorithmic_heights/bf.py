# Bellman-Ford Algorithm

from .helpers import parse_graph
from math import inf, isinf


def nedges(graph):
    count = 0
    for _, outs in graph.items():
        for _ in outs:
            count += 1
    return count


def nodes(graph):
    s = list(graph.keys())
    e = [y["n"] for v in graph.values() for y in v]
    return set(s) | set(e)


def bf(graph, start=1):
    d = {n: inf for n in nodes(graph)}
    d[start] = 0
    for _ in range(nedges(graph) - 1):
        for u, x in graph.items():
            for v in x:
                if d[u] + v["w"] < d[v["n"]]:
                    d[v["n"]] = d[u] + v["w"]
    return ["x" if isinf(v) else v for k, v in d.items()]


def main(file):
    print(*bf(parse_graph(open(file), directed=True, weighted=True)))
