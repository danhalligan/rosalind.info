# Wright-Fisher's Expected Behavior

from .helpers import Parser


def main(file):
    l1, l2 = Parser(file).lines()
    n = int(l1)
    s2 = list(map(float, l2.split()))
    print(*[round(x * n, 3) for x in s2])
