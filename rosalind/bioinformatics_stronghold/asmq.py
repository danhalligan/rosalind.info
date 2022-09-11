# Assessing Assembly Quality with N50 and N75

from .helpers import Parser


def asmq(seqs, n):
    """Assessing Assembly Quality with N50 and N75"""
    lens = sorted([len(x) for x in seqs], reverse=True)
    tot = sum(lens)
    cumsum = 0
    for x in lens:
        cumsum += x
        if cumsum / tot > n / 100:
            return x


def main(file):
    seqs = Parser(file).lines()
    print(asmq(seqs, 50), asmq(seqs, 75))
