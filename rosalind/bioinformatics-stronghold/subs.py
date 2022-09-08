# Finding a Motif in DNA

from .helpers import Parser


def find_motif(seq1, seq2):
    size = len(seq2)
    for i in range(len(seq1) - size + 1):
        if seq1[i : (i + size)] == seq2:
            yield i + 1


def main(file):
    s1, s2 = Parser(file).lines()
    print(*list(find_motif(s1, s2)))
