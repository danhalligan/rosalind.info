# Building a Heap


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


def main(file):
    n, arr = open(file).read().splitlines()
    arr = list(map(int, arr.split()))
    heap = []
    for i, x in enumerate(arr):
        heap.append(x)
        heapify(heap, i)
    print(*heap)
