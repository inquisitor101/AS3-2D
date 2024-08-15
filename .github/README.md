# AS3-2D
This project, code-named **AS3-2D**, stands for: A Simple Structured Solver (AS3) in Two Dimensions (2D). Its purpose is solely for collaborating and learning. As a matter of fact, this endeavor can be traced back to this reddit [post](https://www.reddit.com/r/CFD/comments/1ehllip/collaborating_on_a_cfd_code_for_fun). Also, for active discussion, join the [Slack channel](https://join.slack.com/t/as3-2d/shared_invite/zt-2nxm0hq2u-vwV9I8wIru1YlkN9sUnhQA). 

> [!NOTE]
> As of writing this document, AS3-2D is still under development. Additional information will be provided as the project progresses.

# Overview

AS3-2D uses a nodal discontinuous Galerkin finite element method (DG-FEM) as its spatial discretization strategy. It is restricted to quadrilateral-type elements, utilizing a tensor-product formulation of their 1D counterparts. The temporal discretization uses explicit time-marching schemes only. The solver is capable of handling curvilinear grids with multizone or multiblock configurations.

# Planned Objectives

1. Have an inviscid flow solver, which solves the non-linear Euler equations (EEs).
2. Propose a shared memory parallelization strategy (OpenMP) and implement it. 
3. Propose a distributed-shared hybrid parallelization strategy (MPI+OpenMP) and implement it.
4. Propose a CPU-GPU hybrid parallelization strategy (CUDA-aware MPI + OpenMP) and implement it.
5. Propose explicit SIMD vectorization intrinsics (e.g. SS2, AVX2, AVX512) and implement them.
6. Create a tool that converts a linear Plot3D grid from GMSH, to a high-order grid in native `.as3` format.
7. Utilize the multizone structure to incorporate different buffer/absorbing layers (e.g. PML).
8. Extend the solver to address the compressible Navier-Stokes (NS) equations.

# Getting Started

* To contribute, please read the following [document](./contributing.md).
* For an overview of the coding style used in AS3-2D, see [guidelines](./coding_style.md).
* Documentation will be provided as the project matures.


