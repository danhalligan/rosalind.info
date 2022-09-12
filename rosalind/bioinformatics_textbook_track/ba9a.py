# Construct a Trie from a Collection of Patterns

from collections import defaultdict


def edge_match(graph, node, label):
    m = [x["n"] for x in graph[node] if x["l"] == label]
    return m[0] if m else None


def trie(seqs):
    graph = defaultdict(list)
    count = 1
    for seq in seqs:
        curr_node = 0
        for s in seq:
            match = edge_match(graph, curr_node, s)
            if match:
                curr_node = match
            else:
                graph[curr_node].append({"n": count, "l": s})
                curr_node = count
                count += 1
    return graph


def main(file):
    seqs = open(file).read().splitlines()
    g = trie(seqs)
    for n1, v in g.items():
        for n2 in v:
            print(f"{n1}->{n2['n']}:{n2['l']}")
