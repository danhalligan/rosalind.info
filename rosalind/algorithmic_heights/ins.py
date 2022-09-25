# Insertion Sort

# As mentioned in the problem, we will implement insertion sort and then
# count swaps.

from .helpers import ints


def insertion_sort(x):
    swaps = 0
    for i in range(1, len(x)):
        k = i
        while k > 0 and x[k] < x[k - 1]:
            x[k - 1], x[k] = x[k], x[k - 1]
            swaps += 1
            k -= 1
    return swaps


def main(file):
    _, x = open(file).read().splitlines()
    print(insertion_sort(ints(x)))
