# Enumerating Oriented Gene Orderings

from .helpers import Parser
from itertools import permutations, product


def sign(n):
    """Enumerating Oriented Gene Orderings"""
    perm = list(permutations(range(1, n + 1)))
    sign = list(product([-1, 1], repeat=n))
    res = []
    for p in perm:
        for s in sign:
            res.append([x * y for x, y in zip(s, p)])
    return res


def main(file):
    n = Parser(file).ints()[0]
    res = sign(n)
    print(len(res))
    for i in res:
        print(*i)
