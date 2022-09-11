# Introduction to Set Operations

from builtins import eval
from .helpers import Parser


def main(file):
    n, s1, s2 = Parser(file).lines()
    n = int(n)
    s1, s2 = eval(s1), eval(s2)
    s3 = set(range(1, n + 1))
    print(*[s1 | s2, s1 & s2, s1 - s2, s2 - s1, s3 - s1, s3 - s2], sep="\n")
