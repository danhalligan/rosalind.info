# Distances in Trees

import re
from collections import defaultdict
from math import inf
from heapq import heappop, heappush
from .helpers import nodes


# we will push and pop from a stack of nodes as we work through the
# newick tree backwards (backwards because we can read node names before we
# need to create them!).
#
# "(" corresponds to ending a group of nodes (when working backwards)
# - create an edge from parent to node (move out)
# - move current node to parent and parent to parent's parent
#
# ")" corresponds to starting a group of nodes:
# - set parent to be current node (move in)
# - read node name or take next available integer
# - create a new node with no children
#
# A non-comma token corresponds to a node name
# If it doesn't directly follow a ")" it is a terminal node
# - append to current nodes children
#
# We will create an undirected graph (connecting children to parents)
# and remove the root when we're done (to facilitate a search for shortest
# path...)
def parse_newick(newick, directed=True):
    newick = re.sub(",,", ",.,", newick)
    newick = re.sub(r"\(,", "(.,", newick)
    newick = re.sub(r",\)", ",.)", newick)
    newick = re.sub(r"\(\)", "(.)", newick)
    newick = re.sub(r"^\((.+)\);", r"\1", newick)
    m = re.finditer(r"(\(|[A-z_.]+|,|\))", newick)
    tokens = [x.group() for x in m]

    count = 0
    node_stack = ["0"]
    g = defaultdict(list)
    i = len(tokens) - 1
    while i >= 0:
        if tokens[i] == "(":
            node_stack = node_stack[:-1]
        elif tokens[i] == ")":
            if i + 1 < len(tokens) and tokens[i + 1] not in ",)":
                node = tokens[i + 1]
            else:
                count += 1
                node = str(count)
            g[node_stack[-1]].append({"n": node, "w": 1})
            if not directed:
                g[node].append({"n": node_stack[-1], "w": 1})
            node_stack += [node]
        elif tokens[i] != "," and (i == 0 or tokens[i - 1] != ")"):
            if tokens[i] == ".":
                count += 1
                tokens[i] = str(count)
            g[node_stack[-1]].append({"n": tokens[i], "w": 1})
            if not directed:
                g[tokens[i]].append({"n": node_stack[-1], "w": 1})
        i -= 1
    return g


def dij(graph, start=1):
    dist = {n: inf for n in nodes(graph)}
    dist[start] = 0
    q = []
    heappush(q, (0, start))
    processed = set()

    while q:
        u = heappop(q)[1]
        processed.add(u)
        for v in graph[u]:
            if v["n"] not in processed:
                dist[v["n"]] = min(dist[u] + v["w"], dist[v["n"]])
                heappush(q, (dist[v["n"]], v["n"]))

    return dist


def newick_dist(tree, nodes):
    return dij(parse_newick(tree, directed=False), nodes[0])[nodes[1]]


def main(file):
    contents = open(file).read().split("\n\n")
    if contents[-1] == "":
        contents = contents[:-1]
    trees = [x.split("\n") for x in contents]
    print(*[newick_dist(tree[0], tree[1].split()) for tree in trees])
