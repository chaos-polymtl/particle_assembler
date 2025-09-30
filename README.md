# Particle Assembler

## Overview

Particle Assembler is a software tool designed to generate particle assemblies in cylindrical and cylinder shell geometries using the FIRE (Fast Inertial Relaxation Engine) iterative algorithm. This software is particularly useful for creating initial configurations of particle systems for simulations in computational physics, materials science, and granular mechanics.

## Table of Contents

- [Introduction](#introduction)
- [The FIRE Algorithm](#the-fire-algorithm)
  - [Algorithm Overview](#algorithm-overview)
  - [Mathematical Formulation](#mathematical-formulation)
  - [Algorithm Steps](#algorithm-steps)
- [Geometries](#geometries)
- [Installation](#installation)
- [Usage](#usage)
- [References](#references)
- [License](#license)

## Introduction

Creating well-packed particle assemblies is a fundamental challenge in particle-based simulations. This software implements the FIRE algorithm to efficiently relax particle systems into stable configurations within cylindrical and cylindrical shell geometries.

## The FIRE Algorithm

### Algorithm Overview

The Fast Inertial Relaxation Engine (FIRE) is a molecular dynamics-based minimization algorithm that efficiently finds local energy minima of particle systems. Unlike traditional steepest descent or conjugate gradient methods, FIRE adapts the timestep and introduces artificial damping to accelerate convergence.

The key advantage of FIRE is its ability to handle systems with many degrees of freedom and complex energy landscapes, making it ideal for particle packing problems.

### Mathematical Formulation

The FIRE algorithm is based on a modified molecular dynamics approach with adaptive timestep control. The core equations governing particle motion are:

**Position update:**

$$\mathbf{r}_i(t + \Delta t) = \mathbf{r}_i(t) + \mathbf{v}_i(t) \Delta t + \frac{1}{2}\mathbf{a}_i(t) \Delta t^2$$

**Velocity update:**

$$\mathbf{v}_i(t + \Delta t) = \mathbf{v}_i(t) + \frac{1}{2}[\mathbf{a}_i(t) + \mathbf{a}_i(t + \Delta t)] \Delta t$$

where $\mathbf{r}_i$ is the position of particle $i$, $\mathbf{v}_i$ is the velocity, and $\mathbf{a}_i$ is the acceleration computed from forces:

$$\mathbf{a}_i = \frac{\mathbf{F}_i}{m_i}$$

**Modified velocity (FIRE-specific):**

At each timestep, velocities are modified using:

$$\mathbf{v}_i = (1 - \alpha)\mathbf{v}_i + \alpha \hat{\mathbf{F}}_i |\mathbf{v}_i|$$

where:
- $\alpha$ is a mixing parameter (typically initialized to 0.1)
- $\hat{\mathbf{F}}_i = \mathbf{F}_i / |\mathbf{F}_i|$ is the normalized force direction

**Power criterion:**

The algorithm monitors the power:

$$P = \sum_i \mathbf{F}_i \cdot \mathbf{v}_i$$

When $P > 0$, the system is moving downhill in energy, and the algorithm increases the timestep and decreases $\alpha$ to accelerate convergence. When $P \leq 0$, the system is moving uphill, indicating it has overshot, so velocities are set to zero and parameters are reset.

### Algorithm Steps

The FIRE algorithm proceeds as follows:

1. **Initialize** parameters:
   - $\Delta t = \Delta t_{\text{init}}$ (initial timestep)
   - $\alpha = \alpha_{\text{start}}$ (typically 0.1)
   - $N_{\text{min}} = 5$ (minimum number of steps before increasing $\Delta t$)

2. **MD integration step**:
   - Compute forces $\mathbf{F}_i$ on all particles
   - Update positions and velocities using velocity-Verlet integration
   - Apply velocity modification: $\mathbf{v}_i = (1 - \alpha)\mathbf{v}_i + \alpha \hat{\mathbf{F}}_i |\mathbf{v}_i|$

3. **Check power** $P = \sum_i \mathbf{F}_i \cdot \mathbf{v}_i$:
   
   **If $P > 0$:**
   - Increment counter: $N_{\text{stable}} = N_{\text{stable}} + 1$
   - If $N_{\text{stable}} > N_{\text{min}}$:
     - Increase timestep: $\Delta t = \min(f_{\text{inc}} \Delta t, \Delta t_{\text{max}})$ (with $f_{\text{inc}} = 1.1$)
     - Decrease mixing: $\alpha = f_{\text{dec}} \alpha$ (with $f_{\text{dec}} = 0.99$)
   
   **If $P \leq 0$:**
   - Set all velocities to zero: $\mathbf{v}_i = 0$
   - Decrease timestep: $\Delta t = f_{\text{dec}} \Delta t$ (with $f_{\text{dec}} = 0.5$)
   - Reset mixing: $\alpha = \alpha_{\text{start}}$
   - Reset counter: $N_{\text{stable}} = 0$

4. **Convergence check**:
   - Continue until maximum force $\max_i |\mathbf{F}_i|$ falls below a tolerance threshold
   - Typical convergence criterion: $\max_i |\mathbf{F}_i| < F_{\text{tol}}$

5. **Repeat** steps 2-4 until convergence

## Geometries

### Cylindrical Geometry

Particles are confined within a cylinder of radius $R$ and height $H$. The confinement is enforced through repulsive forces when particles approach the boundaries:

$$\mathbf{F}_{\text{boundary}} = k (d - d_0) \hat{\mathbf{n}}$$

where $d$ is the distance from the boundary, $d_0$ is the threshold distance, $k$ is the stiffness, and $\hat{\mathbf{n}}$ is the normal direction.

**Radial constraint:** For a particle at position $(x, y, z)$:

$$r = \sqrt{x^2 + y^2} \leq R$$

**Axial constraints:**

$$0 \leq z \leq H$$

### Cylinder Shell Geometry

For a cylindrical shell (annular cylinder), particles are confined between an inner radius $R_{\text{inner}}$ and outer radius $R_{\text{outer}}$:

$$R_{\text{inner}} \leq \sqrt{x^2 + y^2} \leq R_{\text{outer}}$$

with the same axial constraint:

$$0 \leq z \leq H$$

## Installation

```bash
# Clone the repository
git clone https://github.com/chaos-polymtl/particle_assembler.git
cd particle_assembler

# Build instructions will be added when source code is available
```

## Usage

```bash
# Usage instructions will be added when source code is available
```

## References

1. Bitzek, E., Koskinen, P., Gähler, F., Moseler, M., & Gumbsch, P. (2006). *Structural Relaxation Made Simple*. Physical Review Letters, 97(17), 170201. [DOI: 10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201)

2. Guénolé, J., Nöhring, W. G., Vaid, A., Houllé, F., Xie, Z., Prakash, A., & Bitzek, E. (2020). *Assessment and optimization of the fast inertial relaxation engine (fire) for energy minimization in atomistic simulations and its implementation in lammps*. Computational Materials Science, 175, 109584.

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.
