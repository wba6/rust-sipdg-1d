# 1D SIPDG Solver in Rust (Sturm–Liouville)

A lightweight, DG-specific, **1D Symmetric Interior Penalty Discontinuous Galerkin (SIPDG)** solver written in **Rust**, targeting **Sturm–Liouville-type elliptic problems** of the form

-(p(x)u')' + q(x)u = f(x)

with Poisson \(-u'' = f\) as a special case. The project is designed to be **transparent and reproducible**: it shows the full computational pipeline from weak formulation → discretization → operator construction → boundary conditions → iterative solve → verification/convergence studies.

This repo is the code component of an Honors project (Spring 2026) focused on clean architecture, DG-oriented data structures, and safe parallelism in Rust.

---

## Goals

- Provide an **educational, method-specific** implementation (not a large general FEM framework) that keeps the DG workflow visible and understandable.
- Support **low-order DG elements** (e.g., linear/quadratic) with **SIPDG fluxes + penalty** for stability and symmetry.
- Maintain solver properties (self-adjoint/coercive assumptions) so the discrete system is **symmetric positive definite**, enabling reliable iterative solvers.
- Offer **matrix-free operator application** (element-by-element + interface-by-interface) and explore **lightweight preconditioning**.
- Use Rust’s safety model + Rayon-style parallel loops to explore **safe parallel DG computation**.

---

## Problem class

- 1D interval domain
- Sturm–Liouville elliptic operator (Poisson as a special case)
- Boundary conditions drawn from **Robin-type** (Dirichlet/Neumann as special cases)

