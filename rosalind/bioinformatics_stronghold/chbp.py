# Character-Based Phylogeny


def pick_informative(characters):
    for i, ch in enumerate(characters):
        if sum(ch) == 2:
            return i, 1
        if sum(ch) == len(ch) - 2:
            return i, 0


def drop_uninformative(characters):
    return [x for x in characters if 1 < sum(x) < len(x) - 1]


def flatten(x):
    if isinstance(x, (list, tuple)):
        return "(" + ",".join(flatten(e) for e in x) + ")"
    else:
        return str(x)


def chbp(names, characters):
    while characters:
        i, v = pick_informative(characters)
        ind = [i for i, x in enumerate(characters[i]) if x == v]
        names[ind[0]] = tuple(names[i] for i in ind)
        del names[ind[1]]
        for ch in characters:
            del ch[ind[1]]
        characters = drop_uninformative(characters)
    return tuple(names)


def main(file):
    names, *characters = open(file).read().splitlines()
    names = names.split()
    characters = [[int(x) for x in list(ch)] for ch in characters]
    res = chbp(names, characters)
    print(flatten(res), ";", sep="")
