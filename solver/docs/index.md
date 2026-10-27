# AS3-2D

**A Simple Structured Solver in Two Dimensions**

AS3-2D is a collaborative computational fluid dynamics project for learning and
experimentation. It solves the compressible Euler equations on structured
quadrilateral grids using a nodal discontinuous Galerkin method and explicit
time integration. Curved elements and multizone grids are supported.

## Start here

- \ref getting_started "Getting started" — build the solver and run an example.
- \ref diagrams "Class hierarchy and dependencies" — explore inheritance, class relationships, and header dependencies.
- \ref architecture "Solver architecture" — follow the main components from grid import to time advancement.
- \ref face_connectivity "Face connectivity" — understand faces, batches, and flattened indices.

## Explore the code

- [Class list](annotated.html) — browse the documented classes and structures.
- [Inheritance diagrams](inherits.html) — see base classes and their implementations.
- [Class hierarchy](hierarchy.html) — browse inheritance as a tree.
- [Source files](files.html) — read the documented headers and implementations.

Individual class pages include inheritance diagrams where applicable. Collaboration
diagrams show relationships such as containment and references between documented
classes. SVG diagrams support zooming and panning.

## Numerical building blocks

| Area | Where to start |
| --- | --- |
| Simulation setup and execution | CDriver, CConfig |
| Multizone grid geometry | CMultizoneGeometry, CSinglezoneGeometry |
| Euler residuals and zone solvers | CMultizoneSolver, ISolver, CEESolver |
| Time integration | ITemporal, CSSPRK3Temporal, CIteration |
| Boundary conditions and interfaces | IBoundary, IInterface |
| Tensor operations | ITensorProduct, CTensorProduct |
| Initial conditions | IInitialCondition |
| Visualization output | IFileVTK, CLegacyBinaryVTK |

## Project

AS3-2D is under active development. The documentation follows the `main` branch
and is rebuilt when changes are pushed to it.

[Source repository](https://github.com/inquisitor101/AS3-2D) ·
[Report an issue](https://github.com/inquisitor101/AS3-2D/issues)
