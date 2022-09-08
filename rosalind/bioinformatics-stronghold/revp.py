# Locating Restriction Sites

from .helpers import Parser


def reverse_pallindromes(seq):
    comp = seq.translate(str.maketrans("ACGT", "TGCA"))
    n = len(seq)
    for i in range(n):
        for size in range(2, 7):
            lo = i - size
            hi = i + size
            if lo < 0 or hi > n or seq[lo:hi] != comp[lo:hi][::-1]:
                break
            else:
                yield [i - size + 1, size * 2]


def main(file):
    seq = Parser(file).fastas()[0].seq
    res = sorted(reverse_pallindromes(seq))
    for row in res:
        print(*row)
