![http://rosalind.info](/images/logo.png)

This repository contains solutions to bioinformatics coding challenges from
[rosalind.info]. Problems are organised by the various different
[locations]:

* [Python Village]: initial problems to learn a few basics about the Python
  programming language.
* [Bioinformatics Stronghold]: problems to discover the algorithms underlying a
  variety of bioinformatics topics.
* [Bioinformatics Armory]: unlike the stronghold in the Armory we solve
* [Bioinformatics Textbook Track]: problems associated with [Bioinformatics
  Algorithms: An Active Learning Approach].
* [Algorithmic Heights]: exercises to accompany the book [Algorithms].
  problems using existing tools.

## Running the solutions

Solutions for each problem are located in individual files inside the directory
for each location.

## Versioning

Versioning of dependencies is controlled here with [poetry]. You can install
the versions of dependencies used here with:

```{shell}
poetry install
```

To run solutions within this environment run, e.g.:

```{shell}
poetry run python3 -m python-village.ini2 ~/Downloads/rosalind_ini2.txt
```

[poetry]: https://python-poetry.org/

## Testing

`pytest-snapshot` is used to test solutions to problems. In many cases solutions
generated will and should exactly match the "Sample Output" given at
rosalind.info. In cases cases, where e.g. ordering is not important, the
expected solutions (in `tests/expected`) have been updated to match code used
here, but are equally valid solutions.

To run the tests use:

```{shell}
poetry run pytest rosalind
```

To update the tests (adding or modifying snapshots / expected output) use:

```{shell}
poetry run pytest --snapshot-update
```

Note that some solutions (that use Entrez) require an email address. This
should be set as an environment variable, e.g.:

```{shell}
export ENTREZ_EMAIL=rosalind.franklin@cam.ac.uk
```

## About

* My rosalind profile: https://rosalind.info/users/danhalligan/

## Locations

### Python Village

The "Python Village" is described as follows.

> If you are completely new to programming, try these initial problems to learn
> a few basics about the Python programming language. You'll get familiar with
> the operations needed to start solving bioinformatics challenges in the
> Stronghold.

<details>
<summary>View progress & solutions</summary>

- [x] [INI1: Installing Python](python-village/ini1.py)
- [x] [INI2: Variables and Some Arithmetic](python-village/ini2.py)
- [x] [INI3: Strings and Lists](python-village/ini3.py)
- [x] [INI4: Conditions and Loops](python-village/ini4.py)
- [x] [INI5: Working with Files](python-village/ini5.py)
- [x] [INI6: Dictionaries](python-village/ini6.py)
</details>

[rosalind.info]: https://rosalind.info
[locations]: https://rosalind.info/problems/locations/
[Python Village]: https://rosalind.info/problems/list-view/?location=python-village
[Bioinformatics Stronghold]: https://rosalind.info/problems/list-view/
[Algorithmic Heights]: https://rosalind.info/problems/list-view/?location=algorithmic-heights
[Bioinformatics Armory]: https://rosalind.info/problems/list-view/?location=bioinformatics-armory
[Bioinformatics Textbook Track]: https://rosalind.info/problems/list-view/?location=bioinformatics-textbook-track
[Bioinformatics Algorithms: An Active Learning Approach]: https://www.bioinformaticsalgorithms.org/
[Algorithms]: https://www.google.co.uk/books/edition/Algorithms/DJSUCgAAQBAJ?hl=en
