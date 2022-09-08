# Find the Longest Path in a DAG

from math import inf
from copy import deepcopy
from collections import defaultdict


# A "reverse" graph, from each node to all incoming nodes
def incoming(graph):
    x = defaultdict(list)
    k = list(graph.keys())
    k += [x["n"] for v in graph.values() for x in v]
    for g in list(set(k)):
        for node in graph[g]:
            x[node["n"]] += [{"n": g, "w": node["w"]}]
    return x


# Nodes with no incoming edges.
def sources(graph):
    inc = incoming(graph)
    return [g for g in graph if len(inc[g]) == 0]


def topological_order(graph):
    graph = deepcopy(graph)
    order = []
    candidates = sources(graph)
    while candidates:
        n = candidates[0]
        order.append(n)
        candidates.remove(n)
        graph[n] = []
        candidates = list(set(sources(graph)) - set(order))
    return order


def longest_path(graph, src, sink):
    score = {}
    path = defaultdict(list)
    for node in graph:
        score[node] = -inf
    score[src] = 0
    path[src] = [src]
    inc = incoming(graph)
    order = topological_order(graph)
    for node in order[order.index(src) + 1 :]:
        if len(inc[node]):
            scores = [score[pre["n"]] + pre["w"] for pre in inc[node]]
            score[node] = max(scores)
            i = scores.index(max(scores))
            path[node] = path[inc[node][i]["n"]] + [node]
    return score[sink], path[sink]


def parse_graph(graph):
    g = defaultdict(list)
    for x in graph:
        n, o = x.split("->")
        node, weight = o.split(":")
        g[n] += [{"n": node, "w": int(weight)}]
    return g


def main(file):
    source, sink, *graph = open(file).read().splitlines()
    score, path = longest_path(parse_graph(graph), source, sink)
    print(score)
    print(*path, sep="->")
