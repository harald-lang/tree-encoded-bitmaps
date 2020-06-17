# Tree-Encoded Bitmaps

This repository contains the source code of
[*Tree-Encoded Bitmaps*](http://db.in.tum.de/~lang/papers/tebs.pdf).

## Quick start

### Check out the code and run the tests
```
git clone https://gitlab.db.in.tum.de/hl/tree-encoded-bitmaps.git teb
git clone https://github.com/peterboncz/bloomfilter-bsd.git dtl
git clone https://gitlab.db.in.tum.de/hl/fastbit.git fastbit
cd teb
mkdir build-debug
cd build-debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j 16 tester
./tester
cd ..
```

### Compression experiments with real-world data sets.
```
mkdir build-release
cd build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 16 real
./real
```

### Experiments with synthetic data sets.

Further experiments, using synthetic data, can be found in the `experiments/` directory.
Each experiment has its own main function and build target, prefixed with `ex_`.
Before an experiment can be executed, the corresponding binary needs to be executed once with the
 environment variable `GEN_DATA=1` set.
This will populate a SQLite database with randomly generated bitmaps, required by the experiment.
All (synthetic) experiments print their CSV results to stdout, which is supposed to be redirect
 into a file.
For details on the output format, please refer to `experiments/{compression,performance}/common.hpp`.

```
make -j 16 ex_compression_uniform
GEN_DATA=1 ./ex_compression_uniform > ex_compression_uniform.out 
```
