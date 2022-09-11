# Compute the Edit Distance Between Two Strings


def edit(s1, s2):
    m = {}
    for j in range(len(s2) + 1):
        m[j, 0] = j
    for i in range(len(s1) + 1):
        m[0, i] = i

    for j in range(len(s2)):
        for i in range(len(s1)):
            if s1[i] == s2[j]:
                m[j + 1, i + 1] = m[j, i]
            else:
                m[j + 1, i + 1] = min([m[j + 1, i], m[j, i], m[j, i + 1]]) + 1

    return m[len(s2), len(s1)]


def main(file):
    s1, s2 = open(file).read().splitlines()
    print(edit(s1, s2))
