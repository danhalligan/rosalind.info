# Adapt SmallParsimony to Unrooted Trees

from collections import defaultdict
from .ba7f import ba6f, root, print_edges


def root_tree(graph, node):
    for child in graph[node]:
        if node in graph[child]:
            graph[child].remove(node)
        root_tree(graph, child)


# parse input (and convert to a rooted tree, setting root to first non-leaf node)
def parse_input(handle):
    n = next(handle)
    n = int(n)
    seqs = {}
    graph = defaultdict(list)
    for i, edge in enumerate(range(n)):
        next(handle)
        f, t = next(handle).rstrip().split("->")
        graph[int(f)].append(i)
        seqs[i] = t

    lines = handle.readlines()
    root = int(lines[0].rstrip().split("->")[0])
    for edge in lines:
        f, t = edge.rstrip().split("->")
        graph[int(f)].append(int(t))

    root_tree(graph, root)
    return seqs, graph


def main(file):
    seqs, graph = parse_input(open(file))
    total_score, seqs = ba6f(graph, seqs)
    print(total_score)
    print_edges(graph, seqs, root(graph))
