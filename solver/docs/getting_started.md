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

## Binary visualization output

Set `OUTPUT_VIS_FORMAT = VTK_XML_BINARY` to write `.vtu` snapshots for ParaView
or VTK. The writer appends uncompressed binary arrays in the machine's byte
order. It caches the fixed mesh in memory and reuses point-data buffers across
snapshots. Each snapshot contains the mesh and the fields selected by
`WRITE_VIS_VARIABLE`. `VTK_LEGACY_BINARY` still selects the `.vtk` writer.

Both writers default to `OUTPUT_VIS_PRECISION = SINGLE` (32-bit coordinates and
point fields). Set `OUTPUT_VIS_PRECISION = DOUBLE` to write 64-bit coordinates
and point fields. This setting controls the file representation. The solver
continues to compute in its usual precision.

Both formats store the physical simulation time in seconds as a double-precision
`TimeValue` field in each snapshot. Open `vis/flow.pvd` for XML output or
`vis/flow.vtk.series` for legacy output in ParaView to animate at the recorded
times. The index is updated after each snapshot and uses paths relative to its
directory, so it can be moved together with the snapshots. See the
[ParaView file-series documentation](https://docs.paraview.org/en/latest/UsersGuide/dataIngestion.html#handling-temporal-file-series).

The XML writer uses 32-bit point indices and cell offsets with 64-bit byte
offsets for appended arrays. Use a VTK-aware reader for the raw binary
payload. See the [VTK XML format specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html).

To compare both writers using ParaView's Python, run from the repository root:

```sh
pvpython --no-mpi solver/tests/vtk_output_test.py solver/bin/release-omp1-ltooff-native0/AS3
```

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
