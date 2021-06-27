# Introduction
An implementation of [Scalable Parallel Graph Coloring](https://www.sciencedirect.com/science/article/pii/S074373150700144X?casa_token=mPUWd9pAQXEAAAAA:SdFYJ60Vx0fmz3cXWIM1TQ1oFWpNpwb2KTAAEvvzzbbum1ybLI-HdrUUKj61Ab0GpXbRhdiVeA) algorithm.

# Requirements
- Open-MPI
- METIS

# Build
Execute this from the project root, be sure you have Open-MPI and METIS installed.
```shell
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
```

Execute this from the project root, there are some sample graphs
provided in directory `graphs`.
# Run
``` shell
# Run with 4 processes using default superstep size (100)
mpirun -np 4 build/pgc graph/rand_1k_5.gr

# Run sequentially using default superstep size
mpirun -np 1 build/pgc graph/rand_1k_5.gr

# Run with 4 processes with superstep size 150
mpirun -np 4 build/pgc graph/rand_1k_5.gr 150
```
