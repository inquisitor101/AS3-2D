# Getting started {#getting_started}

## Requirements

The solver requires a C++20 compiler, GNU Make, and LP64 LAPACK/BLAS libraries.
Parallel builds also require OpenMP. The build options and compiler paths are
configured in `solver/Makefile.gnu`, `solver/Makefile.clang`, and
`solver/Makefile.intel`.

## Build on Linux with GCC

From the repository root:

```sh
cd solver
make -f Makefile.gnu -j8
```

The default build uses OpenMP and writes the executable to
`solver/bin/gnu/release-omp1-ltooff-native0/AS3`.

For a serial build, use `OPENMP=0`. Each build configuration has a separate
output directory; a serial executable has `omp0` in its path.

## Build on macOS with Homebrew LLVM

The Clang makefile defaults to Homebrew LLVM and libomp under `/opt/homebrew/opt`.
Install LLVM, libomp, and LP64 LAPACK/BLAS libraries before building; adjust
`LLVM_PREFIX`, `OMP_PREFIX`, or `LDLIBS` if your installation uses different paths.

```sh
cd solver
make -f Makefile.clang -j8
```

The default executable is
`solver/bin/release-omp1-ltooff-native0/AS3`.

## Run the curved-grid example

Run from the example directory so the mesh paths in `param.cfg` resolve correctly.
Starting from the repository root, for the default macOS Clang build:

```sh
cd example/curved_p5
mkdir -p vis sol
OMP_NUM_THREADS=4 ../../solver/bin/release-omp1-ltooff-native0/AS3 param.cfg
```

For the default GCC build, use
`../../solver/bin/gnu/release-omp1-ltooff-native0/AS3` instead.

The configuration selects the grid, polynomial order, initial condition, boundary
conditions, time integration, and output settings. Visualization output uses the
`vis/flow` prefix in this example. Edit `param.cfg` to change the run settings.

## Build this documentation locally

Install Doxygen and Graphviz, then run from the repository root:

```sh
mkdir -p build/docs
doxygen Doxyfile
python3 .github/doxygen/build_graph_gallery.py build/docs/html
```

Open `build/docs/html/index.html`. The generated HTML is ignored by Git;
only the documentation sources, configuration, and theme files are tracked.

Continue with \ref architecture "Solver architecture" or the
[inheritance diagrams](inherits.html).
