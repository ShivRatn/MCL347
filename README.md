# Course Project: Finite Element Analysis (FEA) for Heat Conduction

This repository contains a set of MATLAB scripts developed as a course project for a Finite Element Analysis (FEA) or Heat Transfer course (likely **MCL347**). The project focuses on solving the steady-state heat conduction problem (Laplace equation) in three different two-dimensional and axisymmetric geometries using the **Finite Element Method (FEM)** with linear triangular elements.

***

## 1. Project Overview

The project implements the FEM from scratch for three distinct geometries, each with specific boundary conditions and coordinate systems, and includes a comparison with known analytical solutions for validation.

| File Name | Geometry | Coordinate System | Boundary Conditions | Key Feature |
| :--- | :--- | :--- | :--- | :--- |
| `rectangle.m` | Rectangular Domain (2D) | Cartesian ($x, y$) | Dirichlet BCs on three sides (one sinusoidal) and implied Neumann on one side. | Heat Flux (Gradient) Calculation |
| `annular_ring.m` | Annular Ring (Axisymmetric 2D) | Polar ($r, \theta$) | Angular-dependent Dirichlet BC on the outer radius. | Axisymmetric FEA Formulation |
| `cylinder.m` | Cylindrical Surface (2D $\theta-z$ plane) | Cylindrical ($\theta, z$) | Z-dependent and angular-dependent Dirichlet BCs on top/bottom faces. | Periodic Boundary Conditions Enforcement |

***

## 2. Mathematical Model and Implementation

### Governing Equation

All simulations solve the steady-state heat conduction equation, which simplifies to the **Laplace equation** in regions with no internal heat generation:

$$\nabla^2 T = 0$$

### Finite Element Formulation

The weak form of the equation is solved using **3-node linear triangular elements (T3)**.

* The **Stiffness Matrix** ($K$) for each element is derived from the integral of the gradient terms ($\nabla N_i \cdot \nabla N_j$) over the element area. The element force vector ($F$) is set to zero as there is no volumetric heat source.
* **Isoparametric mapping** is used to calculate the shape function derivatives and the element area.

---

### `annular_ring.m` - Axisymmetric Heat Transfer 

The code models heat transfer in an annular plate ($r_1 \le r \le r_2$) using a structured polar mesh.

* **Axisymmetric Stiffness:** The standard 2D stiffness matrix ( $\mathbf{K}_{el}$ ) is scaled by the mean radius to incorporate the volume integration in a polar system:
    $$\mathbf{K}_{\text{el, axisymmetric}} = 2\pi r_{\text{avg}} \mathbf{K}_{el}$$

* **Boundary Conditions:**
    * Inner Radius ($r_1 = 0.5$ m): Constant temperature $T_1 = 100^\circ$C.
    * Outer Radius ($r_2 = 1.0$ m): Angular-dependent temperature $T(\theta) = 200 + 50 \sin(\theta)$.

* **Verification:** The numerical solution is verified against the analytical solution for the 2D Laplace equation in polar coordinates, which consists of a constant term and an angular-dependent term.

---

### `rectangle.m` - 2D Cartesian Heat Conduction 🌡️

This module models heat transfer in a **$1.0 \times 1.0$ m square domain** using a structured Cartesian mesh.

* **Boundary Conditions:** Three edges are held at fixed temperatures (Dirichlet BCs). The top edge has a sinusoidal temperature distribution, $T(x) = T_1 (1 + \sin(\pi x / W))$. The fourth side has an implied **Neumann (insulated or zero heat flux) BC**.
* **Post-Processing:** The `calculateGradients` function is used to calculate the **temperature gradients ($\frac{\partial T}{\partial x}, \frac{\partial T}{\partial y}$) which represent the heat flux ($\mathbf{q} \propto -\nabla T$)** at each node, averaging the contribution from adjacent elements.
* **Verification:** The FEM solution is compared against the analytical solution derived using the separation of variables method for the rectangular domain.

---

### `cylinder.m` - Cylindrical Surface Heat Transfer ($\theta-z$ Plane) 🌀

This simulation solves the 2D heat equation in the $\theta-z$ plane, which is relevant for a thin cylinder or when modeling temperature variations on the surface of a solid cylinder of radius $R$.

* **Stiffness Matrix:** The cylindrical stiffness term $\frac{1}{R^2} \frac{\partial^2 T}{\partial \theta^2}$ is explicitly included in the element matrix formulation.
* **Periodic Boundary Conditions:** The code explicitly enforces periodic boundary conditions (BCs) by *merging* the degrees of freedom (DoFs) at $\theta=0$ and $\theta=2\pi$. This ensures temperature continuity across the wrap-around seam of the cylinder.
* **Boundary Conditions:** Angular-dependent Dirichlet BCs are applied on the top ($z=L$) and bottom ($z=-L$) faces: $T(z=\pm L, \theta) = T_1 \pm T_0 \cos(m\theta)$.
* **Verification:** The FEM results are compared to the analytical solution for a solid cylinder with these boundary conditions, involving hyperbolic sine/cosine functions.

***

## 3. Getting Started

### Prerequisites

* MATLAB (R2018a or newer recommended)

### Running the Simulations

The primary function in each file initiates the mesh generation, assembly, solving, and plotting steps.

1. Open the desired file in MATLAB.
2. Run the primary function from the command window:

| File Name | Command to Run |
| :--- | :--- |
| `annular_ring.m` | `axisymmetric_heat_transfer()` |
| `rectangle.m` | `fem_heat_conduction()` |
| `cylinder.m` | `fem_cylindrical_laplace()` |

### Customization

The input parameters can be adjusted at the beginning of each primary function:

* **`annular_ring.m`:** Adjust radii (`r1`, `r2`), boundary temperatures (`T1`), and mesh resolution (`nr`, `ntheta`).
* **`rectangle.m`:** Adjust dimensions (`W`, `H`), base temperature (`T1`), and mesh resolution (`nx`, `ny`).
* **`cylinder.m`:** Adjust cylinder dimensions (`R`, `L`), temperatures (`T1`, `T0`), angular mode number (`m`), and mesh resolution (`nr`, `ntheta`, `nz`).

***

## 4. Expected Output

Each simulation produces multiple figures for visualization and validation:

1. **FEM Solution Plot:** A contour plot of the calculated temperature field.
2. **Analytical Solution Plot:** A contour/surface plot of the exact analytical solution for comparison.
3. **Boundary Condition Verification Plot:** (In `annular_ring.m`) A plot comparing the FEM result to the expected BC values on the outer boundary.
4. **Heat Flux Plot:** (In `rectangle.m`) The contour plot includes superimposed **quiver plots** showing the direction and relative magnitude of the heat flux vectors.
5. **3D Surface Plot:** (In `cylinder.m`) A 3D view of the temperature field mapped onto the cylindrical surface.
6. **2D $\theta-z$ Contour Plot:** (In `cylinder.m`) A 2D contour plot of the temperature field in the $\theta-z$ domain, explicitly showing isotherms.
