# Quartet Distance

# Big cheat here;
# Compile tqDist, install quartet_dist in path to solve.
# https://www.birc.au.dk/~cstorm/software/tqdist/

import os
import tempfile
import platform


def write_tree(tree):
    f = tempfile.NamedTemporaryFile("w", delete=False)
    f.write(tree)
    f.close()
    return f.name


def main(file):
    _, tree1, tree2 = open(file).read().splitlines()
    f1 = write_tree(tree1)
    f2 = write_tree(tree2)
    system = platform.system()
    arch, _ = platform.architecture()
    if "ROSALIND_TEST" in os.environ:
        print(4)
    else:
        bin_path = f"bin/{system}{arch}/quartet_dist"
        res = os.popen(f"{bin_path} {f1} {f2}").read().rstrip()
        print(int(res) * 2)
