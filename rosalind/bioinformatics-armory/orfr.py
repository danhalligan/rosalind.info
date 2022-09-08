# Finding Genes with ORFs

from Bio.Seq import Seq
import re


def trim_seq(x):
    return x[: (len(x) // 3 * 3)]


def translations(x):
    for i in range(3):
        x = trim_seq(x[i:])
        yield str(x.translate())
        yield str(x.reverse_complement().translate())


def main(file):
    best = ""
    seq = str(open(file).read().rstrip())
    for x in translations(Seq(seq)):
        matches = re.findall(r"M[^\*]+", x)
        best = max(matches + [best], key=len)
    print(best)
