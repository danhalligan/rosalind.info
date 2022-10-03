# Identifying Reversing Substitutions

# Some things to note:
# - might not be immediate, e.g. A -> T -> T -> A
# - they can overlap, e.g. A -> T -> A -> T
# - if they occur on internal edges, they will be shared

import re
from .helpers import Parser
from .nwck import parse_newick
from .alph import simplify_tree, nodes


def reverse_graph(graph):
    rev = {}
    for node in graph.keys():
        for child in graph[node]:
            rev[child] = node
    return rev


# The approach here will be to do a DFS. At each position, we find all
# descendants and construct a line of our character table.
def ancestors(tree):
    rev = reverse_graph(tree)
    res = {}
    for x in rev.keys():
        res[x] = []
        node = x
        while node != "0":
            node = rev[node]
            res[x] += [node]
    return res


def extract_position(graph, seqs, pos):
    return {leaf: seqs[leaf][pos] for leaf in nodes(graph) - set("0")}


# Matches any character then not that character (with a negative lookahead) Then
# we match the next character (.) and allow it to repeat (\2*) before matching
# the first again (\1) at the end of the string ($).
def recurrent(node, anc, seq):
    m = re.search(r"(.)(?!\1)((.)\3*)\1$", seq)
    if m:
        g = m.groups()
        return anc[node][len(g[1]) - 1], node, f"{g[0]}->{g[2]}->{g[0]}"


def paths(tree, pos):
    def recurse(node, seq):
        yield recurrent(node, anc, seq + pos[node])
        if node in tree:
            for child in tree[node]:
                yield from recurse(child, seq + pos[node])

    anc = ancestors(tree)
    return recurse("0", "")


def rsub(tree, seqs, i):
    pos = extract_position(tree, seqs, i)
    pos["0"] = ""
    for path in paths(tree, pos):
        if path:
            print(path[0], path[1], i + 1, path[2])
    return None


def main(file):
    handle = open(file)
    tree = parse_newick(next(handle))
    tree = simplify_tree(tree)
    seqs = Parser(handle).fastas()
    seqs = {x.id: x.seq for x in seqs}
    n = len(list(seqs.values())[0])
    for i in range(n):
        rsub(tree, seqs, i)
