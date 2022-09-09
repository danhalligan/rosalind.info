# Enumerating Gene Orders

from itertools import permutations
from .helpers import Parser


def main(file):
    n = Parser(file).ints()[0]
    perm = list(permutations(range(1, n + 1)))
    print(len(perm))
    for i in perm:
        print(*i)
