# Counting Phylogenetic Ancestors

from .helpers import Parser


def main(file):
    print(Parser(file).ints()[0] - 2)
