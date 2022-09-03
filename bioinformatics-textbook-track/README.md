# Bioinformatics textbook track solutions

Solutions to bioinformatic coding challenges that feature within
[Bioinformatics Algorithms: An Active Learning Approach].

The solutions here are relatively naive (and often do not incorporate the
algorithmic improvements suggested within the "Charging station" sections).

The list of problems are also [available online] at [rosalind.info].

[available online]: https://rosalind.info/problems/list-view/?location=bioinformatics-textbook-track
[Bioinformatics Algorithms: An Active Learning Approach]: https://www.bioinformaticsalgorithms.org/
[rosalind.info]: https://rosalind.info

## Running the solutions

Each problem is solved in a separate file. To run this tool to solve a
particular problem, you can either use the `main` function within that file,
or alternatively run `solve.py` from the command line, providing a suitable
input file. For example:

```
./solve.py rosalind_ba1a.txt
```

To reproduce the results for the test data provided at the rosalind website,
you can use the test files included within this repository. For example:

```
./solve.py test/rosalind_ba1a.txt
```
