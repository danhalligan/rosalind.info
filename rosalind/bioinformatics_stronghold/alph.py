# Alignment-Based Phylogeny

# We can solve this with the "Small Parsimony" algorithm

from .helpers import Parser
from .nwck import parse_newick
from math import inf


def nodes(graph):
    s = list(graph.keys())
    e = [y for v in graph.values() for y in v]
    return set(s) | set(e)


# return all leaves of a simple graph
def leaves(graph):
    return nodes(graph) - set(graph.keys())


def extract_position(graph, seqs, pos):
    chars = {}
    for n in nodes(graph) - leaves(graph):
        chars[n] = ""
    for leaf in leaves(graph):
        chars[leaf] = seqs[leaf][pos]
    return chars


def traceback(skp, node, ind):
    bases = ["A", "C", "T", "G", "-"]
    chars = {}
    chars[node] = bases[ind]
    for k, v in skp[node][ind].items():
        if k in skp:
            chars = chars | traceback(skp, k, v)
    return chars


def small_parsimony(graph, chars):
    bases = ["A", "C", "T", "G", "-"]
    sk = {}  # minimum parsimony score of the subtree over possible labels
    skp = {}  # pointer to selected base for each child over possible labels
    to_process = nodes(graph)

    # # initialise leaves
    for leaf in leaves(graph):
        sk[leaf] = [0 if chars[leaf] == c else inf for c in bases]
        to_process.remove(leaf)

    # iterate over available nodes till all are processed
    while to_process:
        for n in list(to_process):
            if all(v in sk for v in graph[n]):
                sk[n], skp[n] = [], []
                for k in bases:
                    tot = 0
                    ptr = {}
                    for d, sk_child in [(d, sk[d]) for d in graph[n]]:
                        score = []
                        for i, c in enumerate(bases):
                            score += [sk_child[i] + (0 if c == k else 1)]
                        tot += min(score)
                        ptr[d] = score.index(min(score))
                    skp[n] += [ptr]
                    sk[n] += [tot]
                to_process.remove(n)

    # Recover sequence
    node = "0"
    score = min(sk[node])
    return score, traceback(skp, node, sk[node].index(score))


def alph(tree, seqs, i):
    # initialise sequences
    for n in nodes(tree) - leaves(tree):
        seqs[n] = ""

    n = len(seqs[list(leaves(tree))[0]])
    total_score = 0
    for pos in range(n):
        chars = extract_position(tree, seqs, pos)
        score, tbchars = small_parsimony(tree, chars)
        total_score += score
        for k, v in tbchars.items():
            seqs[k] += v

    return total_score, seqs


def simplify_tree(graph):
    return {k: [x["n"] for x in v] for k, v in graph.items()}


def main(file):
    handle = open(file)
    tree = parse_newick(next(handle))
    tree = simplify_tree(tree)
    seqs = Parser(handle).fastas()
    seqs = {x.id: x.seq for x in seqs}
    total_score, seqs = alph(tree, seqs, 1)
    print(total_score)
    for node in tree.keys():
        if node != "0":
            print(f">{node}")
            print(seqs[node])
