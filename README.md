# Modelling Reaction-Diffusion in Self-Healing Materials Using FEM

**Course:** Finite Element Methods (ES 622)
**Institution:** Indian Institute of Technology Gandhinagar
**Supervisor:** Prof. Sushobhan Sen

**Submitted by: Group G-11**
* Md. Riaz Ali Mamud (24250055)
* Mayank Awasthi (24250054)
* Nikhil Kaser (24310043)

---

## Project Abstract

This project employs the finite element method (FEM) to model the reaction-diffusion process in a self-healing material, specifically focusing on self-healing polymers. The objective is to understand how the diffusivity of the healing agent impacts the healing process within a material domain. The significance of this study lies in its potential application in extending the lifespan of materials used in various engineering fields, such as construction, automotive, and aerospace.

In this work, the reaction-diffusion equation is numerically solved using a custom FEM code implemented in **MATLAB**.

---

## Governing Equation

The reaction-diffusion equation governing the process is:

$$
\dot{C} = D\nabla^{2}C - kC
$$

This can also be expressed as:

$$
\frac{\partial C}{\partial t} = \nabla \cdot (D \nabla C) - kC
$$

Where:
* **$C$**: Concentration of the healing material ($mol/m^3$)
* **$D$**: Diffusion coefficient ($m^2/s$)
* **$k$**: Reaction rate constant ($1/s$)
* **$\nabla^{2}C$**: The Laplacian, capturing spatial diffusion
* **$-kC$**: The reaction component (first-order decay)

---

## Methodology and Scope

* **Discretization:** The governing equation is discretized in space using an iso-parametric formulation (bilinear quadrilateral elements) and in time using an implicit time-stepping scheme (backward Euler method).
* **Convergence:** A mesh independence study is performed to ensure the solution's accuracy and reliability.
* **Stability:** The stability of the time-stepping scheme (the $\theta$-method) is analyzed, comparing Forward Euler ($\theta=0$) and Crank-Nicolson ($\theta=0.5$).
* **Validation:** The FEM code's predictions are validated against:
    1.  Analytical solutions for diffusion problems (from J. N. Reddy's textbook).
    2.  Commercial FEM software (ABAQUS) for 2D steady-state and transient problems with various geometries (rectangular and hexagonal).
* **Case Studies:** Several test cases are explored to analyze the impact of various parameters on healing efficiency, including:
    * Effect of domain length (simulating crack length)
    * Transient volume fraction of the healing agent
    * Effect of temperature on diffusivity (using the Arrhenius relation)
    * Spatially varying diffusivity
