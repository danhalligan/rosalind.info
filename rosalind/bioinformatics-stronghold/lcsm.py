# Finding a Shared Motif

from .helpers import Parser


def lcsm(seqs):
    seqs = sorted(seqs, key=len)

    maxsub = ""
    s0 = seqs.pop()
    s = len(s0)

    for i in range(s):
        for j in range(s - i + 1):
            if j > len(maxsub) and all(s0[i : i + j] in x for x in seqs):
                maxsub = s0[i : i + j]

    return maxsub


def main(file):
    print(lcsm(Parser(file).seqs()))
