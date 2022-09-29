# Identifying Maximal Repeats

from collections import defaultdict
from .helpers import Parser
from .suff import suff


# This is a bit messy, but the idea here is that we take a suffix tree and
# for each node with children, we get the full path (concatenated edge
# labels) from root. This sequence must be a repeat (since the node has
# children).
#
# Out list `edges` in `main()` contains these paths (from each node with
# children) and the (total) number of descendants that node has.
# Now for each possible number of descendants (which is related to the number
# of repeats), we ensure that we have the maximal repeat by excluding any
# sequences that are fully contained inside others with the same number of
# children.


def as_edges(graph):
    for k in sorted(graph):
        for v in graph[k]:
            yield k, v["n"], v["l"]


def count_descendants(T):
    def count(node):
        if node not in T:
            d[node] = 0
            return 0
        if node not in d:
            d[node] = len(T[node]) + sum(count(x["n"]) for x in T[node])
            return d[node]

    d = {}
    count(0)
    return d


def reverse_graph(graph):
    rev = {}
    for node in graph.keys():
        for child in graph[node]:
            rev[child["n"]] = {"n": node, "l": child["l"]}
    return rev


def pathtoroot(rev, node):
    path = ""
    while node in rev:
        path = rev[node]["l"] + path
        node = rev[node]["n"]
    return path


def internal_edges(graph, node):
    rev = reverse_graph(graph)
    desc = count_descendants(graph)
    if desc[node] > 0:
        yield desc[node], pathtoroot(rev, node)
        for child in graph[node]:
            yield from internal_edges(graph, child["n"])


def as_graph(suff):
    def build_graph(suff, T, n1):
        n2 = n1
        for edge in sorted(suff):
            n2 += 1
            T[n1].append({"n": n2, "l": edge})
            n2 = build_graph(suff[edge], T, n2)
        return n2

    T = defaultdict(list)
    build_graph(suff, T, 0)
    return T


def main(file):
    seq = Parser(file).line()
    tree = as_graph(suff(seq + "$"))
    edges = list(internal_edges(tree, 0))
    d = defaultdict(list)
    for n, edge in edges:
        if len(edge) >= 20:
            d[n].append(edge)
    for repeats in d.values():
        maximal = [x for x in repeats if sum(x in y for y in repeats) == 1]
        print(*maximal, sep="\n")
