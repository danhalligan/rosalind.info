# Sorting by Reversals

from .helpers import Parser
from .rear import sort


def main(file):
    s, t = Parser(file).lines()
    s = list(map(int, s.split()))
    t = list(map(int, t.split()))
    nr, c = sort(s, t)
    print(nr)
    for r in c[0]["p"]:
        print(*r)
