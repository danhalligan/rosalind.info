import yaml
from importlib import resources


def memoize(f):
    cache = {}

    def wrapper(*args):
        if args not in cache:
            cache[args] = f(*args)
        return cache[args]

    return wrapper


def resource_file(file):
    return resources.files("rosalind.resources").joinpath(file)


def yaml_resource(file):
    with open(resource_file(file)) as stream:
        return yaml.safe_load(stream)


@memoize
def genetic_code():
    return yaml_resource("genetic_code.yaml")


@memoize
def codons():
    return yaml_resource("codons.yaml")


@memoize
def aa_mass():
    return yaml_resource("aa_mass.yaml")


@memoize
def blosum62():
    lines = resource_file("blosum62.txt").splitlines()
    header = lines[0].split()
    return dict([x[0], dict(zip(header, map(int, x.split()[1:])))] for x in lines[1:])


@memoize
def pam250():
    lines = resource_file("pam250.txt").splitlines()
    header = lines[0].split()
    return dict([x[0], dict(zip(header, map(int, x.split()[1:])))] for x in lines[1:])
