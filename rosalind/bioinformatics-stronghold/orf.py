# Open Reading Frames

import re
from .helpers import Parser, Dna


def orf(seq):
    for x in [seq, seq.revc()]:
        for i in range(3):
            subseq = x[i : len(x) - (len(x) - i) % 3]
            for m in re.finditer("(?=(M[^\\*]*)\\*)", subseq.translate().seq):
                yield m.group(1)


def main(file):
    seq = Dna(Parser(file).fastas()[0].seq)
    orfs = list(orf(seq))
    orfs = sorted(list(dict.fromkeys(orfs)))
    print("\n".join(orfs))
