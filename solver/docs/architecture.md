# Solver architecture {#architecture}

## Simulation lifecycle

CDriver owns the main simulation components. It reads CConfig, initializes grid
geometry and solver data, advances the solution, and coordinates output.

1. **Read the case.** CConfig parses the input configuration. Grid import supplies
   metadata, coordinates, markers, and interface information.
2. **Construct the geometry.** CMultizoneGeometry owns the zone geometries and
   connectivity. CSinglezoneGeometry describes one zone.
3. **Create the solvers.** CMultizoneSolver owns zone solvers, physical boundary
   conditions, and interfaces. CEESinglezoneSolver implements ISinglezoneSolver for the Euler equations.
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
| ISinglezoneSolver | CEESinglezoneSolver |
| ITemporal | CSSPRK3Temporal |
| IRiemannSolver | CRoeRiemannSolver |
| ITensorProduct | CTensorProduct |
| IBoundary | CSymmetryBoundary, CRiemannStaticMassInletBoundary, CRiemannPressureOutletBoundary |
| IInterface | CEEInterface |
| IInitialCondition | CGaussianPressureIC, CIsentropicVortexIC |
| IFileVTK | CLegacyBinaryVTK, CXMLBinaryVTK |

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

## Local and coupled residual operations

`NResidualComputation::ComputeElementCoupledResidualsCollectiveNowait` traverses
interior, boundary and interface faces. Coupled operations evaluate traces and
exterior states, call the Riemann solver, and lift its numerical flux into the
face-owned buffers. This calculation is independent of the DG formulation.
The existing directional tensor surface operations remain interpolation and
lifting kernels.

After the face-write barrier, `ISinglezoneSolver::ComputeElementLocalResidual`
overwrites one element's local residual. The weak and strong calculations occupy
separate branches. `WEAK_FORM` retains the original volume calculation, including
its interpolation, gradient conversion, weighted contravariant fluxes and
transposed derivative lift. Strong form has separate collocated and non-collocated branches.

For strong form, `CPhysicalElement` precomputes two dense Cartesian
differentiation matrices after the volume metrics are available:
\f[
  (D_x)_{ql} = r_x(q)(D_r)_{ql} + s_x(q)(D_s)_{ql}, \qquad
  (D_y)_{ql} = r_y(q)(D_r)_{ql} + s_y(q)(D_s)_{ql}.
\f]
Rows identify integration/evaluation points and columns identify basis functions.
The reference derivatives are formed from one-dimensional derivative/interpolation
products: dL/dr = dL1D(i,ii)*L1D(j,jj) and
dL/ds = L1D(i,ii)*dL1D(j,jj). The interpolation matrix reduces to the identity
under collocation, up to roundoff. This operator construction also supports
noncollocated nodes, with `nInt2D` rows and `nSol2D` columns. These matrices are
allocated only for strong form and are stored as
`mDerLagrangeInt2Dx` and `mDerLagrangeInt2Dy`.

The collocated strong branch first computes the unweighted Cartesian flux arrays at all
collocated points. The previous dense 2D contraction remains commented in
`CEESinglezoneSolver.cpp`. The active version uses the shared 1D derivative
matrix: r-derivatives follow each row and s-derivatives follow each column,
with only `nInt1D` terms per contraction. Cartesian metrics at the evaluation
point combine these reference derivatives after differentiation, preserving the
Cartesian-flux operator. Its volume contribution is
\f[
  (R_K)_{vq} = -J_q w_q \sum_l
  \left[(D_x)_{ql}(F_x)_{vl} + (D_y)_{ql}(F_y)_{vl}\right].
\f]
The metrics belong to the derivative evaluation point q. Jacobian and quadrature
weights multiply the completed derivative at that same point. Applying these
weights to the input flux samples before differentiation would also
differentiate the weights and yield a different operator.

After volume assembly, the four local faces add
\f$w_k J_{\partial K,k}(n_x F_x + n_y F_y)\f$, using their outward unit normals
and surface Jacobians. The Cartesian fluxes at the collocated boundary nodes
are reused without interpolation or another Euler flux evaluation.

For non-collocated strong form, nodal Cartesian flux coefficients are differentiated
with the dense rectangular matrices. The weighted divergence at quadrature points
is lifted with the transpose of the value interpolation operator:
`R_volume = -L^T W J (Dx Fx + Dy Fy)`. The same nodal flux polynomials are
interpolated to the four faces and their physical normal fluxes are lifted into
the residual. This is a nodal flux-interpolation formulation: for nonlinear Euler
fluxes, its local boundary flux is the trace of the volume flux polynomial, which
generally differs from evaluating the Euler flux of an interpolated state.
The Riemann coupling continues to use state traces. Nonlinear weak/strong
residual equality is not assumed. Non-collocated quadrature does not require a
diagonal mass matrix; the existing dense inverse-mass application is retained.

On curved elements, differentiating interpolated Cartesian fluxes and then
multiplying by metrics generally differs from differentiating interpolated
contravariant fluxes. The dense prototype explicitly uses the former operator;
collocated weak/strong agreement is checked for affine mappings and constant states.
Curved nonlinear residuals are checked against an independent Cartesian
chain-rule reference, with conservation and free-stream tests retained.

`NResidualComputation::ComputeElementResidual` gathers the coupled face buffers
after the local operation and applies the inverse mass matrix exactly once.
Weak form retains its three volume work arrays. Strong form borrows two nodal
Cartesian flux arrays; non-collocated strong form additionally borrows either
a volume-divergence array or two face arrays. Each strong element owns two dense
differentiation matrices of size `nInt2D * nSol2D`. Non-collocated strong form uses
them directly; collocated strong form retains them for reference comparisons
while its active factorized kernel uses the shared 1D derivative matrix.


Mass initialization branches on `CStandardElement::GetIsCollocated()`. For
collocated quadrature, it sets each inverse diagonal entry directly to
`1 / (J[q] * w2D[q])` in a packed vector of `nSol2D` entries and skips the
interpolation products and general matrix inversion. This applies to both weak
and strong forms, including curved elements with varying Jacobians.
Non-collocated elements retain the full mass assembly and inversion path.
The collocated inverse-mass kernel reads consecutive vector entries and scales
each residual entry without summation or scratch storage. Dense non-collocated
storage remains `nSol2D` by `nSol2D`; the tensor factory selects the matching
application kernel.

To explore the implementation, start with CDriver, then CIteration and
CMultizoneSolver. The [file list](files.html) links to the corresponding sources.
