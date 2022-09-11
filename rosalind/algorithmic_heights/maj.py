# Majority Element

# I would naturally solve this with hashing
# Here though, I've implemented Moore's voting algorithm

from .helpers import ints


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


def main(file):
    _, *arrays = open(file).read().splitlines()
    arrays = [ints(arr) for arr in arrays]
    print(*[maj_vote(arr) for arr in arrays])
