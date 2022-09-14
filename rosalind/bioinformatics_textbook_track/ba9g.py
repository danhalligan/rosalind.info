# Construct the Suffix Array of a String


def suffix_array(text):
    seqs = dict((i, text[i:]) for i in range(len(text)))
    return sorted(seqs.keys(), key=lambda x: seqs[x])


def main(file):
    print(*suffix_array(open(file).read().rstrip()), sep=", ")
