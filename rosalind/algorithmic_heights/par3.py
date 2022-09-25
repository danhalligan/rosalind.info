# 3-Way Partition

from .helpers import ints


def par3(A, low=None, high=None):
    if not low:
        low = 0
    if not high:
        high = len(A) - 1
    val = A[low]
    i = low
    while i <= high:
        if A[i] < val:
            A[i], A[low] = A[low], A[i]
            i += 1
            low += 1
        elif A[i] > val:
            A[i], A[high] = A[high], A[i]
            high -= 1
        else:
            i += 1
    return low, high


def main(file):
    handle = open(file)
    next(handle)
    A = ints(next(handle))
    par3(A)
    print(*A)
