# Construct the Suffix Tree of a String

from rosalind.bioinformatics_stronghold.suff import suff, get_edges


def main(file):
    seq = open(file).read().rstrip()
    tree = suff(seq)
    print(*list(get_edges(tree)), sep="\n")
