# Building a Heap

from .helpers import ints


def heapify(heap, i):
    if i == 0:
        return
    parent = (i - 1) // 2
    child = i
    if heap[parent] > heap[child]:
        heapify(heap, parent)
    else:
        heap[parent], heap[child] = heap[child], heap[parent]
        heapify(heap, parent)


def heap(arr):
    heap = []
    for i, x in enumerate(arr):
        heap.append(x)
        heapify(heap, i)
    return heap


def main(file):
    _, arr = open(file).read().splitlines()
    print(*heap(ints(arr)))
