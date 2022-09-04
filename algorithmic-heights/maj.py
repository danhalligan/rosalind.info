# Majority Element

# I would naturally solve this with hashing
# I've implemented this, Moore's voting algorithm and using Counter (cheating!)

from collections import defaultdict, Counter
from .helpers import ints

# import timeit


# Using python Counter
def maj_counter(arr):
    n = len(arr)
    for key, val in Counter(arr).items():
        if val > n / 2:
            return key
    return -1


# Using defaultdict and hashing
# note, in principle, we can stop early, but this involves checking our
# count at each addition...
def maj_hash(arr):
    count = defaultdict(int)
    n = len(arr)
    for x in arr:
        count[x] += 1
    for k, v in count.items():
        if v > n / 2:
            return k
    return -1


def find_candidate(arr):
    maj_index = 0
    count = 1
    for i in range(len(arr)):
        count += 1 if arr[maj_index] == arr[i] else -1
        if count == 0:
            maj_index = i
            count = 1
    return arr[maj_index]


def is_majority(arr, cand):
    count = 0
    for i in range(len(arr)):
        if arr[i] == cand:
            count += 1
    return count > len(arr) / 2


# Using Moore's voting algorithm
def maj_vote(arr):
    cand = find_candidate(arr)
    return cand if is_majority(arr, cand) else -1


def run(arrays, fn):
    return [fn(arr) for arr in arrays]


def main(file):
    info, *arrays = open(file).read().splitlines()
    arrays = [ints(arr) for arr in arrays]

    # Seems that counter is fastest (perhaps no surprise)
    # but its also kind of cheating...
    # print(timeit.timeit(lambda: run(arrays, maj_counter), number = 500))
    # print(timeit.timeit(lambda: run(arrays, maj_hash), number = 500))
    # print(timeit.timeit(lambda: run(arrays, maj_vote), number = 500))

    print(*[maj_vote(arr) for arr in arrays])
