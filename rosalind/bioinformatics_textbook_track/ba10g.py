# Perform a Multiple Sequence Alignment with a Profile HMM

# I'm not 100% sure about the final alignment score calculated here
# but the optimum path appears correct.

from collections import defaultdict
import numpy as np
from .ba10e import state_labels
from .ba10f import pseudocount_profile_hmm
from math import inf, log


def parse_input(handle):
    seq = next(handle).rstrip()
    next(handle)
    θ, σ = map(float, next(handle).rstrip().split())
    next(handle)
    alphabet = next(handle).split()
    next(handle)
    alignment = np.array([list(x) for x in handle.read().splitlines()])
    return seq, θ, σ, alphabet, alignment


def tprob_as_dict(x):
    g = defaultdict(float)
    n = (x.shape[0] - 3) // 3
    lab = state_labels(n)
    for i in range(x.shape[0]):
        for j in range(x.shape[0]):
            g[lab[i], lab[j]] = x[i][j]
    return g


def eprob_as_dict(x, alphabet):
    g = defaultdict(float)
    n = (x.shape[0] - 3) // 3
    lab = state_labels(n)
    for i in range(x.shape[0]):
        for j, c in enumerate(alphabet):
            g[lab[i], c] = x[i][j]
    return g


# Build the weighted HMM graph based on transition matrix
def hmm_graph(tm, n):
    def add_node(a, b):
        g[a].append({"n": b, "w": tm[a, b]})

    g = defaultdict(list)
    for b in ["I0", "M1", "D1"]:
        add_node("S", b)
    for i in range(n):
        a = f"I{i}"
        for b in [a, f"M{i+1}", f"D{i+1}"]:
            add_node(a, b)
    for i in range(1, n):
        for a in [f"M{i}", f"D{i}"]:
            for b in [f"M{i+1}", f"I{i}", f"D{i+1}"]:
                add_node(a, b)
    for a in [f"I{n}", f"M{n}", f"D{n}"]:
        for b in [f"I{n}", "E"]:
            add_node(a, b)

    return g


# Topological order (based on Figure 10.22)
# yields node and column as a tuple
# n: number of valid columns in alignment
# m: sequence length
def topological_order(n, m):
    yield ("S", 0)
    for j in range(n):
        yield (f"D{j+1}", 0)
    for i in range(m):
        yield ("I0", i + 1)
        for j in range(n):
            for c in ["M", "D", "I"]:
                yield (f"{c}{j+1}", i + 1)
    yield ("E", i + 2)


# Function to find previous node and column possibilities based on current node
# and column (see figure 10.21)
def prev_nodes(node, col, n, m):
    if node[0] == "E":
        return [(f"D{n}", m), (f"M{n}", m), (f"I{n}", m)]
    i = int(node[1:])
    if col == 0:
        if i == 1:
            return [("S", 0)]
        else:
            return [(f"D{i-1}", 0)]
    elif node == "I0":
        if col == 1:
            return [("S", 0)]
        else:
            return [("I0", col - 1)]
    elif node == "M1":
        if col == 1:
            return [("S", 0)]
        else:
            return [("I0", col - 1)]
    elif node[0] == "I":
        if col == 1:
            return [(f"D{i}", 0)]
        else:
            return [(f"D{i}", col - 1), (f"M{i}", col - 1), (f"I{i}", col - 1)]
    elif node[0] == "M":
        if col == 1:
            return [(f"D{i-1}", 0)]
        else:
            return [(f"D{i-1}", col - 1), (f"M{i-1}", col - 1), (f"I{i-1}", col - 1)]
    elif node[0] == "D":
        if i == 1:
            return [("I0", col)]
        else:
            return [(f"D{i-1}", col), (f"M{i-1}", col), (f"I{i-1}", col)]
    else:
        print(f"{node} not handled!")


# Convert structure like {n: [{'n': n2, 'w': w}]} to {n: {n2: w}}
def process_graph(graph):
    return {k: {x["n"]: x["w"] for x in v} for k, v in graph.items()}


def main(file):
    seq, θ, σ, alphabet, alignment = parse_input(open(file))
    tprob, eprob = pseudocount_profile_hmm(θ, σ, alphabet, alignment)
    n = (tprob.shape[0] - 3) // 3
    tprob = tprob_as_dict(tprob)
    eprob = eprob_as_dict(eprob, alphabet)

    graph = hmm_graph(tprob, n)
    order = topological_order(n, len(seq))
    graph = process_graph(graph)

    # results will be a dict indexed by tuple (node and column) pointing to score
    # ptr will be indexed by same tuple and point to previous node
    prev = next(order)
    results = {prev: 0}
    ptr = {prev: (None, None)}
    for node, col in order:
        ptr[(node, col)] = 0
        results[(node, col)] = -inf
        for pnode, pcol in prev_nodes(node, col, n, len(seq)):
            if pcol < col and node != "E":
                emission = eprob[node, seq[col - 1]]
            else:
                emission = 1
            prob = log(graph[pnode][node]) + log(emission) + results[(pnode, pcol)]
            if prob > results[(node, col)]:
                results[(node, col)] = prob
                ptr[(node, col)] = (pnode, pcol)

    # traceback
    path = []
    pos = ("E", len(seq) + 1)
    while pos[0]:
        path += [ptr[pos][0]]
        pos = ptr[pos]

    print(*path[::-1][2:])
