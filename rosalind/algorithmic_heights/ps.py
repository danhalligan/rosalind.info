# Partial Sort

from .helpers import ints
from .hea import heap, heapify
from .hs import sift_down, hs

# To start, we insert the first k elements of the input into a max-heap. Then
# for each remaining elements, we add it and heapify to put it into correct
# location. Now we want to remove the biggest element (the first), which we do
# by swapping the root with the last element and sifting down the new root.
#
# Though not specified in the question, I believe the k smallest elements must
# also be sorted...


def ps(arr, k):
    hea = heap(arr[:k])
    for x in arr[k:]:
        hea.append(x)
        heapify(hea, k)
        hea[0] = hea[k]
        hea = hea[:-1]
        sift_down(hea, 0, k - 1)
    return hea


def main(file):
    _, arr, k = open(file).read().splitlines()
    res = ps(ints(arr), int(k))
    print(*hs(res))
