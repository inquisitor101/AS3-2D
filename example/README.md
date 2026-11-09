# Example configurations

Run each configuration from its example directory; grid and output paths are relative to the working directory. Create `vis/` and `sol/` before running. Polynomial degrees and solution nodes come from mesh metadata; `POLY_ORDER` is not a configuration option.

The existing configurations retain their physical parameters, synchronization intervals, and output settings. Spatial defaults and output precision are now written explicitly. Required options have no implicit default.

## Spatial discretization

| Option | Default | Available values / constraints |
| --- | --- | --- |
| `TYPE_DG_FORMULATION` | `WEAK_FORM` | `WEAK_FORM`, `STRONG_FORM`; both support collocated and non-collocated grids. |
| `TYPE_DG_METHOD` | `NODAL_DG` | `NODAL_DG`, `DG_SEM`. |
| `TYPE_QUADRATURE` | `GAUSS_LEGENDRE` | `GAUSS_LEGENDRE`, `LEGENDRE_GAUSS_LOBATTO`, `NEWTON_COTES`. |
| `BONUS_QUADRATURE_FACTOR` | `0` | Non-negative integer `nb`; nodal DG requests quadrature exactness through degree `(2 + nb)*p`. |
| `TYPE_BUFFER_LAYER` | `( NONE )` | `NONE` is the only supported buffer layer. |

`DG_SEM` requires LGL solution nodes, LGL quadrature, and zero bonus. It uses `p + 1` points in each direction and quadrature exactness through degree `2*p - 1`. Collocation requires matching solution and integration nodes, including their locations.

For nodal DG, the requested solution/integration dimensions must be present in the compiled tensor-product factory. `NEWTON_COTES` supports exactness through degree 7, so it cannot satisfy these degree-5 examples, whose nodal DG target starts at degree 10.

The non-collocated strong form differentiates Cartesian flux polynomials interpolated from solution-node fluxes. Increasing quadrature does not change this into evaluating the nonlinear Euler flux at volume quadrature states.

The T106C directory includes three additional configurations with distinct output prefixes:

| Configuration | Method | Formulation | Collocated |
| --- | --- | --- | --- |
| `param_weak.cfg` | `DG_SEM` | `WEAK_FORM` | Yes |
| `param_strong.cfg` | `DG_SEM` | `STRONG_FORM` | Yes |
| `param_noncollocated_strong.cfg` | `NODAL_DG`, LGL quadrature | `STRONG_FORM` | No |

The original `param.cfg` files use weak form. T106C uses DG-SEM; the other original examples use nodal DG with default Gauss-Legendre quadrature.

## Required input and time integration

- `NUMBER_ZONE`: number of zones in the mesh; `TYPE_SOLVER = ( EE )` and `RIEMANN_SOLVER = ( ROE )` select the currently supported equation and Riemann solver. Zone lists follow the configuration parser's zone expansion rules.
- `MESH_FORMAT = AS3`, `INPUT_GRID_FORMAT = BINARY`, and `GRID_DIRECTORY` select the mesh.
- `TEMPORAL_SCHEME = SSPRK3` is the supported time integrator. `START_TIME`, `TIME_STEP`, `MAX_ITER_TIME`, and `CFL_NUMBER` are required. `TIME_STEP` is the synchronization interval; CFL control can take multiple internal steps within it.
- `TYPE_IC` accepts `GAUSSIAN_PRESSURE` or `ISENTROPIC_VORTEX`. Both require `DISTURBANCE_RATIO`, `DISTURBANCE_WIDTH`, `DISTURBANCE_CENTER`, `FREESTREAM_MACH`, and `FREESTREAM_FLOW_ANGLE` (degrees).

For a quick smoke run, copy a configuration and use `TIME_STEP = 1.0e-9`, `MAX_ITER_TIME = 2`, and `WRITE_VIS_FREQ = 2`. Keep the copy in the example directory or use an absolute `GRID_DIRECTORY`. These limits are for a short execution check, not a converged simulation.

## Boundary conditions

Use mesh marker names and supply the parameters below. Multiple marker entries are allowed. Boundary assignments have no implicit default.

| Keyword | Parameters |
| --- | --- |
| `MARKER_BC_PERIODIC` | `( marker1, marker2, translation_x, translation_y )` |
| `MARKER_BC_SYMMETRY` | `( marker )` |
| `MARKER_BC_FARFIELD` | `( marker, density, velocity_x, velocity_y, pressure )` |
| `MARKER_BC_STATIC_MASS_FLOW_INLET_RIEMANN` | `( marker, density, velocity_x, velocity_y )` |
| `MARKER_BC_PRESSURE_OUTLET_RIEMANN` | `( marker, pressure )` |

## Output

`OUTPUT_VIS_FORMAT` is required and accepts `VTK_LEGACY_BINARY` (`.vtk`) or `VTK_XML_BINARY` (`.vtu`). `OUTPUT_VIS_PRECISION` defaults to `SINGLE` (32-bit); `DOUBLE` (64-bit) is also supported for coordinates and point fields in both formats.

`WRITE_VIS_FREQ` (positive), `OUTPUT_VIS_FILENAME`, `OUTPUT_SOL_FILENAME`, and `WRITE_VIS_VARIABLE` are required. Visualization variables are `DENSITY`, `MOMENTUM`, `TOTAL_ENERGY`, `PRESSURE`, `VELOCITY`, `VORTICITY`, `MACH`, `TEMPERATURE`, and `ENTROPY`.
