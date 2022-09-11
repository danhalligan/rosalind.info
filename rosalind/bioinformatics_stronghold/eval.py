# Expected Number of Restriction Sites

from .helpers import Parser


def eval(n, s, x):
    """Expected Number of Restriction Sites"""
    gc = sum([s.count(x) for x in "GC"])
    p = ((1 - x) / 2) ** (len(s) - gc) * (x / 2) ** gc
    return p * (n - len(s) + 1)


def main(file):
    n, s, a = Parser(file).lines()
    n = int(n)
    a = map(float, a.split())
    print(*[round(eval(n, s, x), 3) for x in a])
