# Using the Spectrum Graph to Infer Peptides

from .helpers import Parser
from .full import weight_graph


def sgra(ions):
    """Using the Spectrum Graph to Infer Peptides"""

    def infer_peptide(w, seq):
        for w2, aa in graph[w]:
            yield from infer_peptide(w2, seq + aa)
        yield seq

    graph = weight_graph(ions)
    w = min(ions)
    return max(list(infer_peptide(w, "")), key=len)


def main(file):
    weights = [float(x) for x in Parser(file).lines()]
    print(sgra(weights), sep="\n")
