import yaml
from importlib import resources
from functools import cache


def resource_file(file):
    return resources.files("rosalind.resources").joinpath(file)


def yaml_resource(file):
    with open(resource_file(file)) as stream:
        return yaml.safe_load(stream)


def matrix_resource(file):
    lines = open(resource_file(file)).read().splitlines()
    header = lines[0].split()
    return {x[0]: dict(zip(header, map(int, x.split()[1:]))) for x in lines[1:]}


@cache
def genetic_code():
    return yaml_resource("genetic_code.yaml")


@cache
def codons():
    return yaml_resource("codons.yaml")


@cache
def aa_mass():
    return yaml_resource("aa_mass.yaml")


@cache
def blosum62():
    return matrix_resource("blosum62.txt")


@cache
def pam250():
    return matrix_resource("pam250.txt")
