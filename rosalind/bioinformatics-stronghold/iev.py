# Calculating Expected Offspring

from .helpers import Parser


def iev(v):
    p = [1, 1, 1, 0.75, 0.5, 0]
    return sum([x[0] * x[1] * 2 for x in zip(v, p)])


def main(file):
    print(iev(Parser(file).ints()))
