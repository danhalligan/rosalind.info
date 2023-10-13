# Construct a Suffix Tree from a Suffix Array

# from .ba9g import suffix_array


def main(file):
    txt, sa, lcp = open(file).read().splitlines()
    sa = [int(x) for x in sa.split(", ")]
    lcp = [int(x) for x in lcp.split(", ")]

    # txt = "panamabananas$"
    # sa = suffix_array(txt)
    # print(sa)
    # lcp = [0, 0, 1, 1, 3, 3, 1, 0, 0, 0, 2, 2, 0, 0]

    edges = []
    lcp.append(0)
    levels = [0]
    for i in range(len(txt)):
        suffix = txt[sa[i] :]
        before = lcp[i]
        after = lcp[i + 1]
        diff = after - before
        if diff == 0:
            edges.append(suffix[before:])
        elif diff > 0:
            edges.append(suffix[after:])
            levels.append(after)
        else:
            edges.append(suffix[before:])
            while True:
                if levels[-2] == after:
                    edges.append(suffix[after : levels[-1]])
                    levels.pop()
                    break
                elif levels[-2] < after:
                    edges.append(suffix[after : levels[-1]])
                    levels.pop()
                    levels.append(after)
                    break
                else:
                    edges.append(suffix[levels[-2] : levels[-1]])
                    levels.pop()
                    continue
    print(*sorted(edges), sep="\n")
