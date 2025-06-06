Production-ready MiniSAT
========================

[![Build Status](https://travis-ci.org/master-keying/minisat.svg?branch=master)](https://travis-ci.org/master-keying/minisat)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/master-keying/minisat?svg=true)](https://ci.appveyor.com/project/horenmar/minisat)

Forked off MiniSAT 2.2, this repository aims at providing
a production-ready version of the famous library.

### Key principles:
- Keep the algorithm intact.
- Adapt for the modern C++ 2011 & 2014 standard.
- Support GCC, LLVM and MSVC.
- Support Linux, Windows native and MinGW/MSYS.
- Fix complier warnings (`-Wall -Wextra`).
- Ensure code quality by automated build and tests.
- Clean the std. output and error output.



Known limitations (aka TODO)
----------------------------

1. `libgz` library was removed. Please use `gunzip | minisat` instead.
2. Memory and CPU limits/stats are not supported well on all platforms.

Keen to help? Feel free to fork!



Directory overview
------------------

- `minisat/mtl`     Mini Template Library
- `minisat/utils`   Generic helper code (I/O, Parsing, CPU-time, etc)
- `minisat/core`    A core version of the solver
- `minisat/simp`    An extended solver with simplification capabilities
- `README.md`       This read-me file
- `LICENSE`         Licence files



Building
--------

### a) Quick way
```bash
./build.sh
```


### b) Traditional way

Release version: without assertions, statically linked, etc

```bash
export MROOT=<minisat-dir> # or `setenv` in cshell
cd minisat
cd { core | simp }
make rs
cp minisat_static <install-dir>/minisat
```



Examples
========

Run minisat with same heuristics as version 2.0:

```minisat <cnf-file> -no-luby -rinc=1.5 -phase-saving=0 -rnd-freq=0.02```
