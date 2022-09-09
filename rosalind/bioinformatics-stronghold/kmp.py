# Speeding Up Motif Finding

from .helpers import Parser


def kmp_preprocess(seq):
    """KMP preprocessing algorithm"""
    j = -1
    b = [j]
    for i in range(len(seq)):
        while j >= 0 and seq[i] != seq[j]:
            j = b[j]
        j += 1
        b.append(j)
    return b[1:]


def main(file):
    seq = Parser(file).fastas()[0].seq
    print(*kmp_preprocess(seq))
