# Counting Point Mutations

from .helpers import Parser


def hamm(s1, s2):
    """Calculate Hamming distance"""
    return sum(xi != yi for xi, yi in zip(s1, s2))


def main(file):
    print(hamm(*Parser(file).lines()))
