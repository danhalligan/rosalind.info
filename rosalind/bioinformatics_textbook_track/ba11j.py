# Find a Highest-Scoring Modified Peptide against a Spectrum

from .ba11c import masses
from itertools import accumulate
import networkx as nx

# We need to consider k changes that lead to pv having same length as sv
# This is basically an alignment problem which we can solve by finding the
# longest path in a DAG (or shortest path if we negate the weights).
# We build the DAG from figure 11.20 and then find best path from source to sink

# Since there's negative weights we need to use Bellman-Ford, and there's lots
# of connections, so I'm using networkx which is quick.

# We could also solve this with a dynamic programming approach as for the
# sequence alignment problems. If I have time, I'll try this approach so
# I don't have to use networkx...


def build_dag_nx(cs, sv, k):
    """Build a directed graph with k layers"""
    dag = nx.DiGraph()

    # Add diagonal edges within each t
    for t in range(k + 1):
        for ix, i in enumerate(cs[:-1]):
            d = cs[ix + 1] - cs[ix]
            for j in range(len(sv) - d):
                dag.add_edge((i, j, t), (i + d, j + d, t), weight=-sv[j + d - 1])

    # Add edges between the t layers (non-diagonals)
    # Each corresponds to a modified peptide
    for t in range(k):
        for ix, i in enumerate(cs[:-1]):
            for j in range(len(sv)):
                for m in range(j + 1, len(sv) + 1):
                    if m - j != cs[ix + 1] - i:
                        dag.add_edge(
                            (i, j, t), (cs[ix + 1], m, t + 1), weight=-sv[m - 1]
                        )

    return dag


def recover_string(peptide, path):
    """Given best path, recover the string of amino acids as required by problem"""
    j = 0
    string = ""
    for c, aa in enumerate(peptide):
        m = masses()[aa]
        d = path[c + 1][1] - (j + m)
        string += aa if d == 0 else f"{aa}({d :+})"
        j = path[c + 1][1]
    return string


def main(file):
    peptide, sv, k = open(file).read().splitlines()
    sv = list(map(int, sv.split()))
    k = int(k)
    mass = [masses()[x] for x in peptide]
    cs = [0] + list(accumulate(mass))  # the cumulative sum of masses
    dag = build_dag_nx(cs, sv, k)

    # TODO
    # We should check all t up to k really (since we might get better score
    # with fewer modifications), but I haven't done that yet
    end = (cs[-1], len(sv), k)
    path = nx.shortest_path(dag, (0, 0, 0), end, weight="weight", method="bellman-ford")
    print(recover_string(peptide, path))
