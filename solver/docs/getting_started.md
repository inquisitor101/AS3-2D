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

## Weak and strong DG formulations

`TYPE_DG_FORMULATION` selects `WEAK_FORM` (the default) or `STRONG_FORM`,
independently of `TYPE_DG_METHOD`. Existing configurations without this option
continue to use weak form. The selected formulation is included in the startup
spatial-discretization report.

Strong form supports both collocated and non-collocated quadrature. For the
supplied LGL grids, select collocated strong DG-SEM with:

```text
TYPE_DG_METHOD = DG_SEM
TYPE_QUADRATURE = LEGENDRE_GAUSS_LOBATTO
BONUS_QUADRATURE_FACTOR = 0
TYPE_DG_FORMULATION = STRONG_FORM
```

The collocated strong local operation uses 1D row/column differentiation of the
Cartesian fluxes and combines them with the Cartesian metrics, then adds outward physical fluxes
on all four faces. The coupled operation adds the Riemann numerical flux
contributions. See \ref architecture for the assembly
and metric conventions.

Run the formulation regressions with:

```sh
bash solver/tests/run_dg_formulation.sh
```

These use a debug OpenMP build with address and undefined-behavior sanitizers.
They check Cartesian differentiation matrices, strong residuals against an
independent chain-rule reference, conservation, free-stream preservation,
SSP-RK3 updates and configuration errors, with one and four threads.
Synthetic skew and curved elements cover orders 1 through 9. The imported curved
example has an existing constant-state residual of approximately `1.5e-8` in both
forms; the test bounds this fixture error separately from the tight weak/strong
constant-state comparison and the near-roundoff synthetic-element checks.

For non-collocated strong form, use `TYPE_DG_METHOD = NODAL_DG` with GL or
LGL quadrature and `TYPE_DG_FORMULATION = STRONG_FORM`. The Cartesian flux
polynomials are defined by Euler flux values at solution nodes. Dense rectangular
derivatives evaluate their divergence at quadrature points, followed by weighted
value lifting into the solution basis. Local face terms trace the same flux
polynomials. Coupled Riemann terms still use interpolated solution traces.
Increasing quadrature integrates these nodal flux polynomials more accurately;
it does not evaluate the nonlinear volume flux afresh at the extra quadrature
points, as the existing weak implementation does.

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
