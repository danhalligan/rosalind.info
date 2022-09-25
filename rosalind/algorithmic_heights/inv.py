# Counting Inversions

from .helpers import ints


def mer2(a1, a2):
    a = []
    s = len(a1)
    count, i = 0, 0
    while a1 and a2:
        if a1[0] <= a2[0]:
            i += 1
            a += [a1.pop(0)]
        else:
            count += s - i
            a += [a2.pop(0)]
    return a + a1 + a2, count


def ms2(a):
    if len(a) > 1:
        mid = len(a) // 2
        a1, c1 = ms2(a[:mid])
        a2, c2 = ms2(a[mid:])
        a, c = mer2(a1, a2)
        return a, c1 + c2 + c
    else:
        return a, 0


def main(file):
    _, array = open(file).read().splitlines()
    _, count = ms2(ints(array))
    print(count)
