# Construct a Profile HMM

# from book:
# ignore columns for which the fraction of space symbols is greater than or
# equal to a column removal threshold q

import numpy as np


def parse_input(handle):
    θ = float(next(handle).rstrip())
    next(handle)
    alphabet = next(handle).split()
    next(handle)
    alignment = np.array([list(x) for x in handle.read().splitlines()])
    return θ, alphabet, alignment


# compute column or row in emission probability matrix
def index(i, type):
    if type == "ins":
        return (i + 1) * 3 + 1
    else:
        return {"match": 0, "del": 1}[type] + 3 * i + 2


# There must be an easier way, but this divides rows by rowMax
# ignore zeros in the row, or when the rowsum is zero.
def rownorm(x, inc_zeros=False, min_val=0.0):
    if inc_zeros and sum(x) == 0:
        x[:] = 1
    with np.errstate(divide="ignore", invalid="ignore"):
        new = x / sum(x)
        new[x == 0.0] = min_val
        return new


def normalise(x, inc_zeros=False, min_val=0.0):
    return np.array([rownorm(r, inc_zeros=inc_zeros, min_val=min_val) for r in x])


# This problem is very fussy about formatting.
# I think we need to ensure correct rounding and certainly tab-delimiting
def print_mat(mat, rl, cl):
    print(*cl, sep="\t")
    for i, row in enumerate(mat):
        r = [rl[i]] + [round(x, 3) if x > 0.0 else "0" for x in row]
        print(*r, sep="\t")


def print_tprob(x):
    n = (x.shape[0] - 3) // 3
    print_mat(x, state_labels(n), state_labels(n))


def print_eprob(x, alphabet):
    n = (x.shape[0] - 3) // 3
    print_mat(x, state_labels(n), alphabet)


def state_labels(n):
    x = ["S", "I0"]
    for i in range(1, n + 1):
        x += [f"M{i}", f"D{i}", f"I{i}"]
    x += "E"
    return x


def transition_mat(n):
    x = np.zeros((n * 3 + 3, n * 3 + 3), dtype=float)
    return x


# Initialise a emission probability matrix
def emission_mat(n, m):
    return np.zeros((n * 3 + 3, m), dtype=float)


def profile_hmm(θ, alphabet, alignment):
    valid_col = np.mean(alignment == "-", axis=0) < θ
    valid_len = sum(valid_col)
    end = valid_len * 3 + 2
    tprob = transition_mat(valid_len)
    eprob = emission_mat(valid_len, len(alphabet))

    for seq in alignment:
        pind = 0
        j = -1
        for i, char in enumerate(seq):
            if valid_col[i]:
                j += 1
                if char == "-":
                    ind = index(j, "del")
                else:
                    ind = index(j, "match")
                tprob[pind, ind] += 1
                pind = ind
            else:
                if char != "-":
                    ind = index(j, "ins")
                    tprob[pind, ind] += 1
                    pind = ind
            if char != "-":
                eprob[ind, alphabet.index(char)] += 1
        tprob[pind, end] += 1

    tprob = normalise(tprob)
    eprob = normalise(eprob)

    return tprob, eprob


def main(file):
    θ, alphabet, alignment = parse_input(open(file))
    tprob, eprob = profile_hmm(θ, alphabet, alignment)
    print_tprob(tprob)
    print("--------")
    print_eprob(eprob, alphabet)
