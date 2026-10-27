# Solver architecture {#architecture}

## Simulation lifecycle

CDriver owns the main simulation components. It reads CConfig, initializes grid
geometry and solver data, advances the solution, and coordinates output.

1. **Read the case.** CConfig parses the input configuration. Grid import supplies
   metadata, coordinates, markers, and interface information.
2. **Construct the geometry.** CMultizoneGeometry owns the zone geometries and
   connectivity. CSinglezoneGeometry describes one zone.
3. **Create the solvers.** CMultizoneSolver owns zone solvers, physical boundary
   conditions, and interfaces. CEESolver implements ISolver for the Euler equations.
4. **Advance in time.** ITemporal defines time integration. CSSPRK3Temporal
   orchestrates SSP-RK3 stages; CIteration computes spatial residuals and advances
   elements within a stage.
5. **Write results.** COutput coordinates output, with IFileVTK implementations
   providing visualization files.

## Interfaces and implementations

See \ref diagrams "Class hierarchy and dependencies" for the inheritance overview
and class collaboration diagrams. The [full graphical hierarchy](inherits.html)
is also available.
These interfaces are useful starting points:

| Interface | Implementations |
| --- | --- |
| ISolver | CEESolver |
| ITemporal | CSSPRK3Temporal |
| IRiemannSolver | CRoeRiemannSolver |
| ITensorProduct | CTensorProduct |
| IBoundary | CSymmetryBoundary, CRiemannStaticMassInletBoundary, CRiemannPressureOutletBoundary |
| IInterface | CEEInterface |
| IInitialCondition | CGaussianPressureIC, CIsentropicVortexIC |
| IFileVTK | CLegacyBinaryVTK |

Follow a class link to view its methods, source, and available diagrams.
Inheritance diagrams show base and derived classes. Collaboration diagrams also
include implementation relationships between documented classes, such as member
variables that refer to other classes. They describe static relationships in the
code rather than execution order.

## Grid connectivity

Connectivity is grouped by logical direction and by interior, physical-boundary,
and interface faces. Batches collect faces with shared metadata. Flat indices
allow the iteration code to traverse these collections efficiently.

See \ref face_connectivity "Face connectivity" for the data model, index spaces,
and ownership conventions.

## Parallel execution

COpenMP configures shared-memory execution. CIteration::AdvanceStage creates the
OpenMP team and synchronizes surface writes before element updates. Temporal
updates operate on individual elements.

To explore the implementation, start with CDriver, then CIteration and
CMultizoneSolver. The [file list](files.html) links to the corresponding sources.
