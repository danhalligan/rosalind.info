# Merge Two Sorted Arrays

from .helpers import ints


def mer(a1, a2):
    a = []
    while a1 and a2:
        a.append(a1.pop(0) if a1[0] < a2[0] else a2.pop(0))
    return a + a1 + a2


def main(file):
    _, a1, _, a2 = open(file).read().splitlines()
    print(*mer(ints(a1), ints(a2)))
