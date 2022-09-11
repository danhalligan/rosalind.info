# Finding the Longest Multiple Repeat

from collections import defaultdict
from .helpers import Parser


def build_seq(node, rev, edges):
    seq = ""
    while node in rev:
        seq = edges[node] + seq
        node = rev[node]
    return seq


# traverse graph and find longest seq up to a node with a least k leaves?
def lrep(seq, k, graph):
    # parse data and build structures
    edges = defaultdict(list)
    rev = {}
    heads, tails = set(), set()
    for edge in graph:
        n1, n2, i, n = edge.split()
        tails.add(n2)
        heads.add(n1)
        rev[n2] = n1
        edges[n2] = seq[(int(i) - 1) : (int(i) + int(n) - 1)]

    # count the number of descendents (leaves) from each node
    descendents = defaultdict(int)
    for leaf in tails - heads:
        while leaf in rev:
            leaf = rev[leaf]
            descendents[leaf] += 1

    candidates = [x for x in descendents if descendents[x] >= k]
    seqs = [build_seq(cand, rev, edges) for cand in candidates]
    return max(seqs, key=len)


def main(file):
    seq, k, *g = Parser(file).lines()
    print(lrep(seq, int(k), g))
