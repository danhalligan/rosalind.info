# Median

from .helpers import ints
from .par3 import par3


# s and e are the indices of the start and end of the middle partition
def med(A, k):
    def find_med(A, k, low, high):
        s, e = par3(A, low, high)
        if k < s:
            return find_med(A, k, low, s)
        elif k > e:
            return find_med(A, k, e + 1, high)
        else:
            return A[s]

    return find_med(A, k, 0, len(A) - 1)


def main(file):
    handle = open(file)
    _ = next(handle)
    A = ints(next(handle))
    k = int(next(handle))
    print(med(A, k - 1))
