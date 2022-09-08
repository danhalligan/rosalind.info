# Align Two Strings Using Affine Gap Penalties

from .ba5e import blosum62


def insert_indel(word, i):
    return word[:i] + "-" + word[i:]


def global_affine(v, w, sigma=-11, epsilon=-1):
    score = blosum62()
    m, x, y = {}, {}, {}
    px, pm, py = {}, {}, {}

    x[0, 0], m[0, 0], y[0, 0] = 0, 0, 0
    for i in range(1, len(v) + 1):
        x[i, 0] = sigma + (i - 1) * epsilon
        m[i, 0] = sigma + (i - 1) * epsilon
        y[i, 0] = 10 * sigma
    for j in range(1, len(w) + 1):
        y[0, j] = sigma + (j - 1) * epsilon
        m[0, j] = sigma + (j - 1) * epsilon
        x[0, j] = 10 * sigma

    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s = [x[i - 1, j] + epsilon, m[i - 1, j] + sigma]
            x[i, j] = max(s)
            px[i, j] = s.index(x[i, j])

            s = [y[i, j - 1] + epsilon, m[i, j - 1] + sigma]
            y[i, j] = max(s)
            py[i, j] = s.index(y[i, j])

            s = [x[i, j], m[i - 1, j - 1] + score[v[i - 1]][w[j - 1]], y[i, j]]
            m[i, j] = max(s)
            pm[i, j] = s.index(m[i, j])

    i, j = len(v), len(w)
    a1, a2 = v, w

    scores = [x[i, j], m[i, j], y[i, j]]
    max_score = max(scores)
    s = scores.index(max_score)

    while i * j != 0:
        if s == 0:
            if px[i, j] == 1:
                s = 1
            i -= 1
            a2 = insert_indel(a2, j)
        elif s == 1:
            if pm[i, j] == 1:
                i -= 1
                j -= 1
            else:
                s = pm[i, j]
        else:
            if py[i, j] == 1:
                s = 1
            j -= 1
            a1 = insert_indel(a1, i)

    for _ in range(i):
        a2 = insert_indel(a2, 0)
    for _ in range(j):
        a1 = insert_indel(a1, 0)

    return max_score, a1, a2


def main(file):
    s1, s2 = open(file).read().splitlines()
    print(*global_affine(s1, s2), sep="\n")
