# Find the Longest Substring Shared by Two Strings

from collections import defaultdict
from rosalind.bioinformatics_stronghold.suff import suff
import re


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


def leaves(graph):
    tails = list(n["n"] for nodes in graph.values() for n in nodes)
    heads = set(graph.keys())
    return set(tails) - heads


# reverse a (child points to parent), retaining edge labels
def reverse_graph(graph):
    rev = {}
    for node in graph.keys():
        for child in graph[node]:
            rev[child["n"]] = {"n": node, "l": child["l"]}
    return rev


def purple_edges(tree, colors):
    def _purple_edges(node, seq):
        if colors[node] == "purple":
            anyp = any(colors[n["n"]] == "purple" for n in tree[node])
            if anyp:
                for child in tree[node]:
                    yield from _purple_edges(child["n"], seq + child["l"])
            else:
                yield seq

    yield from _purple_edges(0, "")


# initialise leaves. The final nodes in our graph have no contents, but
# here we'll set its color based on the label from the previous nodes
# edges
def leaf_colors(tree):
    rev = reverse_graph(tree)
    colors = defaultdict(lambda: None)
    for leaf in leaves(tree):
        edge = rev[leaf]["l"]
        m = re.search(r"([$#])", edge)
        colors[leaf] = "red" if m.group() == "$" else "blue"
    return colors


# Add "colours" to a paired suffix tree.
# "red" == sequence 1 (only), ending "$"
# "blue" == sequence 2 (only), ending "#"
# "purple == shared.
def color_tree(tree, colors):
    q = list(tree.keys())
    while q:
        for node in q:
            cols = [colors[n["n"]] for n in tree[node]]
            if all(cols):
                cols = set(cols)
                if cols == set(["red"]):
                    colors[node] = "red"
                elif cols == set(["blue"]):
                    colors[node] = "blue"
                else:
                    colors[node] = "purple"
                q.remove(node)

    return colors


def longest_shared_substring(seq1, seq2):
    tree = suff(seq1 + "$" + seq2 + "#")
    tree = as_graph(tree)
    colors = color_tree(tree, leaf_colors(tree))
    return max(purple_edges(tree, colors), key=lambda x: len(x))


def main(file):
    seq1, seq2 = open(file).read().splitlines()
    print(longest_shared_substring(seq1, seq2))
