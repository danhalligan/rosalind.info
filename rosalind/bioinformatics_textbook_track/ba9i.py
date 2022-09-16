# Construct the Burrows-Wheeler Transform of a String

from .ba9g import suffix_array


def bwt(seq):
    return "".join(seq[i - 1] for i in suffix_array(seq))


def main(file):
    print(bwt(open(file).read().rstrip()))
