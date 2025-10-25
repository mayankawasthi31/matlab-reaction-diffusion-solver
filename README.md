# Modelling Reaction-Diffusion in Self-Healing Materials Using FEM

[cite_start]**Course:** Finite Element Methods (ES 622) [cite: 4]
[cite_start]**Institution:** Indian Institute of Technology Gandhinagar [cite: 13]
[cite_start]**Supervisor:** Prof. Sushobhan Sen [cite: 8]

**Submitted by: Group G-11**
* [cite_start]Md. Riaz Ali Mamud (24250055) [cite: 6]
* [cite_start]Mayank Awasthi (24250054) [cite: 6]
* [cite_start]Nikhil Kaser (24310043) [cite: 7]

---

## Project Abstract

[cite_start]This project employs the finite element method (FEM) to model the reaction-diffusion process in a self-healing material, specifically focusing on self-healing polymers[cite: 25]. [cite_start]The objective is to understand how the diffusivity of the healing agent impacts the healing process within a material domain[cite: 26]. [cite_start]The significance of this study lies in its potential application in extending the lifespan of materials used in various engineering fields, such as construction, automotive, and aerospace[cite: 27].

[cite_start]In this work, the reaction-diffusion equation is numerically solved using a custom FEM code implemented in **MATLAB**[cite: 28].

---

## Governing Equation

[cite_start]The reaction-diffusion equation governing the process is[cite: 28]:

$$
\dot{C} = D\nabla^{2}C - kC
$$

This can also be expressed as:

$$
\frac{\partial C}{\partial t} = \nabla \cdot (D \nabla C) - kC
$$

Where:
* [cite_start]**$C$**: Concentration of the healing material ($mol/m^3$) [cite: 29]
* [cite_start]**$D$**: Diffusion coefficient ($m^2/s$) [cite: 29]
* [cite_start]**$k$**: Reaction rate constant ($1/s$) [cite: 29]
* [cite_start]**$\nabla^{2}C$**: The Laplacian, capturing spatial diffusion [cite: 29]
* [cite_start]**$-kC$**: The reaction component (first-order decay) [cite: 29, 91]

---

## Methodology and Scope

* [cite_start]**Discretization:** The governing equation is discretized in space using an iso-parametric formulation (bilinear quadrilateral elements) and in time using an implicit time-stepping scheme (backward Euler method)[cite: 35, 40, 53].
* [cite_start]**Convergence:** A mesh independence study is performed to ensure the solution's accuracy and reliability[cite: 37, 694].
* [cite_start]**Stability:** The stability of the time-stepping scheme (the $\theta$-method) is analyzed, comparing Forward Euler ($\theta=0$) and Crank-Nicolson ($\theta=0.5$)[cite: 53, 729, 733, 734].
* **Validation:** The FEM code's predictions are validated against:
    1.  [cite_start]Analytical solutions for diffusion problems (from J. N. Reddy's textbook)[cite: 38, 698, 699].
    2.  [cite_start]Commercial FEM software (ABAQUS) for 2D steady-state and transient problems with various geometries (rectangular and hexagonal)[cite: 39, 55, 1009, 1469].
* [cite_start]**Case Studies:** Several test cases are explored to analyze the impact of various parameters on healing efficiency[cite: 42, 50, 56, 1613], including:
    * [cite_start]Effect of domain length (simulating crack length) [cite: 43, 1616]
    * [cite_start]Transient volume fraction of the healing agent [cite: 44, 1695]
    * [cite_start]Effect of temperature on diffusivity (using the Arrhenius relation) [cite: 46, 1764, 1771]
    * [cite_start]Spatially varying diffusivity [cite: 47, 1828]
