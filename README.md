![http://rosalind.info](/images/logo.png)

This repository contains solutions to bioinformatics coding challenges from
[rosalind.info]. Problems are organised by the various different
[locations]:

* [Algorithmic Heights]
* [Bioinformatics Armory]
* [Bioinformatics Textbook Track]: problems associated with [Bioinformatics
  Algorithms: An Active Learning Approach].

[rosalind.info]: https://rosalind.info
[locations]: https://rosalind.info/problems/locations/
[Algorithmic Heights]: https://rosalind.info/problems/list-view/?location=algorithmic-heights
[Bioinformatics Armory]: https://rosalind.info/problems/list-view/?location=bioinformatics-armory
[Bioinformatics Textbook Track]: https://rosalind.info/problems/list-view/?location=bioinformatics-textbook-track
[Bioinformatics Algorithms: An Active Learning Approach]: https://www.bioinformaticsalgorithms.org/

## Running the solutions

Solutions for each problem are located in individual files inside the directory
for each location.

To run this tool to solve a particular problem, you can either run the script
(providing the input file as a command line argument) or use `solve.py` which
will solve based on the input file name.

```{shell}
python3 ini.py test/rosalind_fibo.txt
```

```{shell}
./solve.py test/rosalind_fibo.txt
```

## About

* My rosalind profile: https://rosalind.info/users/danhalligan/
