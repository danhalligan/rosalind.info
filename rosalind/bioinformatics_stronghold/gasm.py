# Genome Assembly Using Reads

from .helpers import Parser
from .dbru import dbru
from .pcov import find_cycle, join_cycle


def kmers(seq, k):
    """Return all kmers of length k from sequence seq"""
    return [seq[i : (i + k)] for i in range(len(seq) - k + 1)]


def extract_chain(graph):
    ch1 = find_cycle(graph)
    graph = dict(filter(lambda x: x[0] not in ch1, graph.items()))
    return ch1, graph


def gasm(seqs):
    """Genome Assembly Using Reads"""
    s0 = seqs[0][:-1]
    for k in range(len(s0) + 1, 0, -1):
        subseqs = [i for x in seqs for i in kmers(x, k)]
        graph = dict(dbru(subseqs))
        try:
            ch1, graph = extract_chain(graph)
            ch2, graph = extract_chain(graph)
            if len(graph) == 0:
                return join_cycle(ch1)
        except KeyError:
            continue


def main(file):
    print(gasm(Parser(file).lines()))
