# Inferring Peptide from Full Spectrum

from math import isclose
from collections import defaultdict
from .helpers import Parser
from .spec import match_mass


# Return a dictionary that let's us look up complementary ion pairs
# e.g. "PRO" and "TEIN" for full peptide "PROTEIN"
# We can find pairs, because their weights approximately sum to that
# of the peptide.
def find_pairs(peptide, ions):
    pairs = {}
    for w in ions:
        for w2 in ions:
            if isclose(w2 + w, peptide):
                pairs[w] = w2
    return pairs


# This builds a graph of possibly adjacent ions.
# ions can be considered adjacent in a graph if the difference in
# their masses is close to the mass of a known amino acid.
def weight_graph(ions):
    graph = defaultdict(list)
    for i in range(len(ions)):
        for j in range(i + 1, len(ions)):
            aa = match_mass(ions[j] - ions[i])
            if aa:
                graph[ions[i]] += [[ions[j], aa]]
    return graph


def full(peptide, ions):
    """Inferring Peptide from Full Spectrum"""

    def infer_peptide(w, seq, rm):
        for w2, aa in graph[w]:
            if w2 in rm:
                continue
            if len(seq) + 1 == target_len:
                yield seq + aa
            else:
                yield from infer_peptide(w2, seq + aa, rm + [w, pairs[w]])

    graph = weight_graph(ions)
    pairs = find_pairs(peptide, ions)
    target_len = int(len(ions) / 2 - 1)
    w = min(ions)
    return list(infer_peptide(w, "", [w, pairs[w]]))


def main(file):
    weights = [float(x) for x in Parser(file).lines()]
    print(*full(weights[0], weights[1:]), sep="\n")
