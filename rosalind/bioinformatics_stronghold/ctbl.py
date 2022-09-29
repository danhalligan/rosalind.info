# Creating a Character Table

import re
from .helpers import nodes
from .nwck import parse_newick


# The approach here will be to do a DFS. At each position, we find all
# descendants and construct a line of our character table.
def descendants(T):
    def recurse(node):
        if node not in T:
            d[node] = []
            return []
        if node not in d:
            children = [x["n"] for x in T[node]]
            d[node] = children + [y for x in T[node] for y in recurse(x["n"])]
            return d[node]

    d = {}
    recurse("0")
    return d


# This is hacky and relies on knowing that the taxa are not just numbers (which
# are used for other internal nodes not of interest)
def ctbl(tree):
    desc = descendants(tree)
    allnodes = sorted(x for x in nodes(tree) if not re.match(r"^\d+$", x))
    for subnodes in desc.values():
        if not len(subnodes) or all([x in subnodes for x in allnodes]):
            continue
        yield "".join(["1" if n in subnodes else "0" for n in allnodes])


def main(file):
    tree = parse_newick(open(file).read().rstrip())
    for code in ctbl(tree):
        print(code)
