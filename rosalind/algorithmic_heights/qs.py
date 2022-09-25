# Quick Sort

from .helpers import ints
from .par3 import par3


# Quick Sort is very similar to our median algorithm: recursive partitioning.
def qs(A):
    def qsr(A, low, high):
        s, e = par3(A, low, high)
        if low < s:
            qsr(A, low, s)
        if high > e:
            qsr(A, e + 1, high)

    qsr(A, 0, len(A) - 1)


def main(file):
    handle = open(file)
    _ = next(handle)
    A = ints(next(handle))
    qs(A)
    print(*A)
