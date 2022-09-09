# Matching a Spectrum to a Protein

from .helpers import Parser
from .conv import conv
from .spec import spectrum


def prsm(s, r):
    """Matching a Spectrum to a Protein"""
    mult = [conv(spectrum(x), r)[1] for x in s]
    return [max(mult), s[mult.index(max(mult))]]


def main(file):
    dat = Parser(file).lines()
    n = int(dat[0])
    res = prsm(dat[1 : 1 + n], [float(x) for x in dat[1 + n :]])
    print(*res, sep="\n")
