# Binary Search

from .helpers import ints


def bins(x, arr):
    low = 0
    high = len(arr) - 1
    while low <= high:
        mid = (low + high) // 2
        if x == arr[mid]:
            return mid + 1
        elif x < arr[mid]:
            high = mid - 1
        else:
            low = mid + 1
    return -1


def main(file):
    _, _, arr, m = open(file).read().splitlines()
    arr = ints(arr)
    m = ints(m)
    print(*[bins(x, arr) for x in m])
