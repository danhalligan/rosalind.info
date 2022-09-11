# Linguistic Complexity of a Genome

from .helpers import Parser
from .suff import suff, get_edges


def ling(seq):
    s = sum(len(edge) for edge in get_edges(suff(seq)))
    m = sum(min(4**k, len(seq) - k + 1) for k in range(1, len(seq) + 1))
    return s / m


def main(file):
    seq = Parser(file).line()
    print(ling(seq))
