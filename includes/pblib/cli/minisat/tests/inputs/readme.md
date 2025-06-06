# Notes about the inputs

The basic set of inputs were taken from the [SATLIB benchmark page]
(http://www.cs.ubc.ca/~hoos/SATLIB/benchm.html). For simplicity, they are
split into two main folders, `SAT` and `UNSAT`, depending on whether the
problem is satisfiable or unsatisfiable.

The inputs were split into 2 groups, integration tests and benchmarks.
The integration tests are such inputs that the current implementation
can solve in a reasonable amount of time, and their list is contained
within the `easy.txt` file. The benchmarks are DIMACS inputs that the
current implementation needs at least a minute to solve, and their list
is contained within the `benchmarks.txt` file.

## Input sources

### DIMACS Benchmark Instances (original source: [DIMACS Benchmark set for SAT](ftp://dimacs.rutgers.edu/pub/challenge/satisfiability/benchmarks/cnf/))
* AIM: Artificially generated Random-3-SAT - 48 instances satisfiable, 24 unsatisfiable
* LRAN: Large Random-3-SAT instances - 3 instances, all satisfiable
* JNH: Random SAT instances with variable length clauses - 16 instances satisfiable, 34 instances unsatisfiable
* DUBOIS: Randomly generated SAT instances - 13 instances, all unsatisfiable
* GCP: Large SAT-encoded Graph Colouring problems - 4 instances, all satisfiable
* PARITY: Instances for problem in learning the parity function - 20 instances, all satisfiable
* II: Instances from a problem in inductive inference - 41 instances, all satisfiable
* HANOI: SAT-encoding of Towers of Hanoi - 2 instances, all satisfiable
* BF: Circuit fault analysis: bridge fault - 4 instances, all unsatisfiable
* SSA: Circuit fault analysis: single-stuck-at fault - 4 instances satisfiable, 4 instances unsatisfiable
* PHOLE: Pigeon hole problem - 5 instances, all unsatisfiable
* PRET: Encoded 2-colouring forced to be unsatisfiable - 8 instances, all unsatisfiable
