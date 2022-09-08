# Solve the Turnpike Problem

# ## PROBLEM
#
# We have a set of distances between positions
# e.g. -10 -8 -7 -6 -5 -4 -3 -3 -2 -2 0 0 0 0 0 2 2 3 3 4 5 6 7 8 10

# We we are trying to infer a set of positions that would be
# consistent with this, e.g. 0 2 4 7 10

# It is easy to infer distances from positions: all_differences([0, 2, 4, 7, 10])
# but harder to do the reverse.

# ## SOLUTION OUTLINE
#
# Note that:
# * We can assume that the first position is at 0.
#
# * Logically, the maximum position must be the maximum of the distances (10)
#
# * We can ignore negative and 0-distances since all positions must be unique
#   Thus, we will consider only the distances [2, 2, 3, 3, 4, 5, 6, 7, 8]
#
# * To find a solution, we can prune distances from this list while adding
#   positions to our candidate list till there are no more distances.
#
# * When adding a candidate position, we must check that all distances between
#   this position and our current list are within the known distances list. This
#   can be accomplished by subtracting Counter objects: xdist - dist.
#
# * To proceed we note that the next largest value above is 7. So either, 7
#   exists as a position, or, since we know we have 10 in our list, we may have
#   1 as a position.
#
# * We consider both possibilities, check if the distances implied by adding
#   this are consistent with the known differences, and if so, recurse.
#
# * In our first iteration, a position of 8 is consistent, since we have
#   differences 8-0 = 8 and 10-8 = 2 in our known list. It could also be 2.
#
#   In our next iteration, we cannot have 7 as a position, since this would
#   imply our difference list would contain 10 - 7 == 3, which it does not...

from collections import Counter
from itertools import product


def all_differences(pos):
    return sorted([a - b for a, b in product(pos, repeat=2)])


def difference(pos, d):
    return Counter(abs(x - d) for x in pos)


def turnpike(dist, pos):
    if not dist:
        yield sorted(pos)
    else:
        for x in [max(dist), max(pos) - max(dist)]:
            xdist = difference(pos, x)
            if not (xdist - dist):
                yield from turnpike(dist - xdist, pos + [x])


def main(file):
    dist = open(file).read().split()
    dist = list(map(int, dist))
    pos = [0, max(dist)]
    dist = Counter(filter(lambda x: x > 0, dist[:-1]))
    # for solution in turnpike(dist, pos):
    #     print(*solution)
    print(*next(turnpike(dist, pos)))
