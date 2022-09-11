# Comparing Spectra with the Spectral Convolution

from collections import Counter
from .helpers import Parser


def conv(s1, s2):
    """Comparing Spectra with the Spectral Convolution"""
    x = sorted([round(i - j, 5) for i in s1 for j in s2])
    return Counter(x).most_common()[0]


def main(file):
    l1, l2 = Parser(file).lines()
    s1 = list(map(float, l1.split()))
    s2 = list(map(float, l2.split()))
    res = conv(s1, s2)
    print(res[1], res[0], sep="\n")
