#!/usr/bin/env python3

import typer
import importlib
import re


def main(file: str):
    """Solve problem based on filename"""
    problem = re.search(r".+_(.+)\.txt", file).groups()[0]
    module = importlib.import_module(problem)
    getattr(module, "main")(file)


if __name__ == "__main__":
    typer.run(main)
