# Implement SmallParsimony

from collections import defaultdict
from math import inf


# return all nodes of a simple graph
def nodes(graph):
    s = list(graph.keys())
    e = [y for v in graph.values() for y in v]
    return set(s) | set(e)


# return all leaves of a simple graph
def leaves(graph):
    return set(y for v in list(graph.values()) for y in v if not graph[y])


# return all root node of a simple graph
def root(graph):
    rev = reverse_graph(graph)
    node = list(nodes(graph))[0]
    while node in rev:
        node = rev[node]
    return node


# reverse a simple graph (child points to parent)
def reverse_graph(graph):
    rev = {}
    for node in graph:
        for child in graph[node]:
            rev[child] = node
    return rev


def parse_input(handle):
    n = next(handle)
    n = int(n)
    seqs = {}
    graph = defaultdict(list)
    for i, edge in enumerate(range(n)):
        f, t = next(handle).rstrip().split("->")
        graph[int(f)].append(i)
        seqs[i] = t
    for edge in handle.readlines():
        f, t = edge.rstrip().split("->")
        graph[int(f)].append(int(t))
    return seqs, graph


# print (bidirectional) edges
def print_edges(graph, seqs, node):
    for child in graph[node]:
        dist = sum(a != b for a, b in zip(seqs[node], seqs[child]))
        print(f"{seqs[node]}->{seqs[child]}:{dist}")
        print(f"{seqs[child]}->{seqs[node]}:{dist}")
        print_edges(graph, seqs, child)


def extract_position(graph, seqs, pos):
    chars = {}
    for n in nodes(graph) - leaves(graph):
        chars[n] = ""
    for leaf in leaves(graph):
        chars[leaf] = seqs[leaf][pos]
    return chars


def traceback(skp, node, ind):
    bases = ["A", "C", "T", "G"]
    chars = {}
    chars[node] = bases[ind]
    for k, v in skp[node][ind].items():
        if k in skp:
            chars = chars | traceback(skp, k, v)
    return chars


def small_parsimony(graph, chars):
    bases = ["A", "C", "T", "G"]
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
    node = root(graph)
    score = min(sk[node])
    return score, traceback(skp, node, sk[node].index(score))


def ba6f(graph, seqs):
    # initialise sequences
    for n in nodes(graph) - leaves(graph):
        seqs[n] = ""

    total_score = 0
    for pos in range(len(seqs[0])):
        chars = extract_position(graph, seqs, pos)
        score, tbchars = small_parsimony(graph, chars)
        total_score += score
        for k, v in tbchars.items():
            seqs[k] += v

    return total_score, seqs


def main(file):
    seqs, graph = parse_input(open(file))
    total_score, seqs = ba6f(graph, seqs)
    print(total_score)
    print_edges(graph, seqs, root(graph))
