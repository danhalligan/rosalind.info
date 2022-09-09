# Ordering Strings of Varying Length Lexicographically

import re
from .helpers import Parser
from itertools import product


def lexv(s, n):
    """Ordering Strings of Varying Length Lexicographically"""
    s = ["_"] + s
    perm = list(product(s, repeat=n))
    perm = ["".join(x) for x in perm]
    perm = [re.sub("_+$", "", x) for x in perm]
    return list(filter(lambda x: "_" not in x, perm[1:]))


def main(file):
    l1, l2 = Parser(file).lines()
    print(*lexv(l1.split(), int(l2)), sep="\n")
