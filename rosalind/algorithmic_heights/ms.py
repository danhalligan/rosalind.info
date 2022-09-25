# Merge Sort

from .mer import mer
from .helpers import ints


def ms(a):
    if len(a) > 1:
        mid = len(a) // 2
        return mer(ms(a[:mid]), ms(a[mid:]))
    else:
        return a


def main(file):
    _, a = open(file).read().splitlines()
    print(*ms(ints(a)))
