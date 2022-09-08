# Compute the 2-Break Distance Between a Pair of Genomes

import re
from collections import defaultdict


def find_component(node, graph):
    q = [node]
    visited = set()
    while q:
        node = q.pop(0)
        visited.add(node)
        for n in graph[node]:
            if n not in visited:
                q += [n]
    return visited


# In the format chosen here, I take a directed edge +x (e.g. "+1")
# and assign the node at the “head” of this edge as -1 and the
# node at the “tail” as +1.
def parse_genome_graph(s):
    g = defaultdict(list)
    for component in re.findall(r"\((.+?)\)", s):
        component = list(map(int, component.split()))
        for i in range(len(component) - 1):
            g[component[i]] += [-component[i + 1]]
            g[-component[i + 1]] += [component[i]]
        g[component[-1]] += [-component[0]]
        g[-component[0]] += [component[-1]]
    return g


# merge the graphs (assumes gp and gp have the same nodes)
def breakpoint_graph(p, q):
    bg = {}
    for node in p.keys():
        bg[node] = p[node] + q[node]
    return bg


def ba6c(genomes):
    bg = breakpoint_graph(*genomes)
    nodes = set(bg.keys())
    n_blocks = len(nodes) // 2
    n_components = 0
    while len(nodes):
        res = find_component(next(iter(nodes)), bg)
        nodes = nodes - res
        n_components += 1

    return n_blocks - n_components


def main(file):
    genomes = [parse_genome_graph(s) for s in open(file).read().splitlines()]
    print(ba6c(genomes))
