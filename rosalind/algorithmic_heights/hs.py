# Heap Sort

from .helpers import ints
from .hea import heap


# sift-down the element at heap[start] to its proper place
def sift_down(heap, start, end):
    root = start
    while root * 2 + 1 <= end:
        left = root * 2 + 1
        right = left + 1
        swap = root
        if heap[swap] < heap[left]:
            swap = left
        if right <= end and heap[swap] < heap[right]:
            swap = right
        if swap != root:
            heap[root], heap[swap] = heap[swap], heap[root]
            root = swap
        else:
            return


def hs(x):
    x = heap(x)
    i = len(x) - 1
    while i > 0:
        x[0], x[i] = x[i], x[0]
        i -= 1
        sift_down(x, 0, i)
    return x


def main(file):
    _, arr = open(file).read().splitlines()
    print(*hs(ints(arr)))
