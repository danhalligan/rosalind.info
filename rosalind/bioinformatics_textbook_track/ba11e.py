# Sequence a Peptide

from collections import defaultdict
from .ba11d import prefixes2peptide
from .ba11c import masses
from .ba7a import nodes
from math import inf


# Bellman-Ford Algorithm (required for the negative weights)
# This version finds longest path in a DAG from source to sink and returns path
def bf(start, graph):
    d = {n: -inf for n in nodes(graph)}
    d[start] = 0
    p = {start: []}
    n = nodes(graph)
    for _ in range(len(n) - 1):
        for u, x in graph.items():
            for v in x:
                if d[u] + v["w"] > d[v["n"]]:
                    d[v["n"]] = d[u] + v["w"]
                    p[v["n"]] = p[u] + [v["n"]]
    return d, p


def build_dag(v, masses):
    v = [0] + v
    dag = defaultdict(list)
    m = {v: k for k, v in masses.items()}
    for i in range(len(v)):
        for j in range(i + 1, len(v)):
            if j - i in m:
                dag[i].append({"n": j, "w": v[j]})
    return dag


def main(file):
    m = masses()
    vector = list(map(int, open(file).read().split()))
    dag = build_dag(vector, m)
    _, path = bf(0, dag)
    print(prefixes2peptide(path[len(vector)], m))
