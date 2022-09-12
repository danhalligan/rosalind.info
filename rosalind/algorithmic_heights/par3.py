# 3-Way Partition


def par3(A):
    start = 0
    end = len(A) - 1
    val = A[0]
    i = 0
    while i <= end:
        if A[i] < val:
            A[i], A[start] = A[start], A[i]
            i += 1
            start += 1
        elif A[i] > val:
            A[i], A[end] = A[end], A[i]
            end -= 1
        else:
            i += 1
    return start, end


def main(file):
    handle = open(file)
    next(handle)
    A = list(map(int, next(handle).split()))
    par3(A)
    print(*A)
