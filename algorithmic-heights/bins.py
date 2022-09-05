# Binary Search
import sys


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
    n, m, arr, ints = open(file).read().splitlines()
    arr = list(map(int, arr.split()))
    ints = list(map(int, ints.split()))
    print(*[bins(x, arr) for x in ints])


if __name__ == "__main__":
    main(sys.argv[1])
