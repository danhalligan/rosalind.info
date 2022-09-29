# Newick Format with Edge Weights

import re
from collections import defaultdict
from .nwck import dij


# A slightly hacky modification of nwck:parse_newick that only works if
# weights are present.
# TODO: rewrite to be a single function...


def parse_newick(newick, directed=True):
    newick = re.sub(",,", ",.,", newick)
    newick = re.sub(r"\(,", "(.,", newick)
    newick = re.sub(r",\)", ",.)", newick)
    newick = re.sub(r"\(\)", "(.)", newick)
    newick = re.sub(r"^\((.+)\);", r"\1", newick)
    m = re.finditer(r"(\(|([A-z_.]*:\d+)|,|\))", newick)
    tokens = [x.groups()[0] for x in m]

    count = 0
    node_stack = ["0"]
    g = defaultdict(list)
    i = len(tokens) - 1
    while i >= 0:
        if tokens[i] == "(":
            node_stack = node_stack[:-1]
        elif tokens[i] == ")":
            if i + 1 < len(tokens) and tokens[i + 1] not in ",)":
                if tokens[i + 1][0] == ":":
                    weight = tokens[i + 1][1:]
                    count += 1
                    node = str(count)
                else:
                    node, weight = tokens[i + 1].split(":")
            g[node_stack[-1]].append({"n": node, "w": int(weight)})
            if not directed:
                g[node].append({"n": node_stack[-1], "w": int(weight)})
            node_stack += [node]
        elif tokens[i] != "," and (i == 0 or tokens[i - 1] != ")"):
            if tokens[i] == ".":
                count += 1
                tokens[i] = str(count)
            node, weight = tokens[i].split(":")
            g[node_stack[-1]].append({"n": node, "w": int(weight)})
            if not directed:
                g[node].append({"n": node_stack[-1], "w": int(weight)})
        i -= 1
    return g


def newick_dist(tree, nodes):
    return dij(parse_newick(tree, directed=False), nodes[0])[nodes[1]]


def main(file):
    contents = open(file).read().split("\n\n")
    if contents[-1] == "":
        contents = contents[:-1]
    trees = [x.split("\n") for x in contents]
    print(*[newick_dist(tree[0], tree[1].split()) for tree in trees])
