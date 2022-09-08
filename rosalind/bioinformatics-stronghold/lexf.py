# Enumerating k-mers Lexicographically

from .helpers import Parser
from itertools import product


def main(file):
    l1, l2 = Parser(file).lines()
    set = l1.split(" ")
    n = int(l2)
    perm = ["".join(x) for x in product(set, repeat=n)]
    print(*sorted(perm), sep="\n")
