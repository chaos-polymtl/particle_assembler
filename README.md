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

### Requirements

- C++ compiler with C++17 support
- CMake 3.10 or higher
- Python 3.x with NumPy and PyVista (for visualization)

### Build Instructions

```bash
# Clone the repository
git clone https://github.com/chaos-polymtl/particle_assembler.git
cd particle_assembler

# Create build directory
mkdir build
cd build

# Configure and build
cmake ..
make

# The executable will be created as particle_assembler
```

## Usage

### Running the Particle Assembler

The particle assembler generates a packed particle assembly in a cylindrical shell geometry:

```bash
# Run the assembler
./particle_assembler
```

This will:
1. Initialize N particles randomly within the cylindrical shell
2. Gradually grow particle radii from an initial low packing fraction to the target packing fraction
3. Use FIRE algorithm to relax overlaps at each growth step
4. Output the final particle configuration to `output.txt`

NOTE: The particle_assembler currently does not have an API to parametrize it from input file. This will be added in the near future.

### Configuration

Key parameters are defined in `source/particle_assembler.cc`:

- `N`: Number of particles (default: 70486)
- `phi_target`: Target packing fraction (default: 0.30)
- `r_in`: Inner radius of cylindrical shell (default: 0.0064 m)
- `thickness`: Shell thickness (default: 0.0238 - r_in m)
- `height`: Cylinder height (default: 0.25 m)
- `k_pair`: Particle-particle interaction stiffness (default: 1e3)
- `k_wall`: Particle-wall interaction stiffness (default: 1e3)
- `grow_rate`: Particle radius growth rate per cycle (default: 1.02)
- `fire_max_steps`: Maximum FIRE iterations per cycle (default: 100000)
- `fire_dt`: FIRE timestep (default: 1e-5)
- `FIRE_FTOL`: Force convergence tolerance (default: 1e-6)

### Output Format

The software outputs a file `output.txt` with the following format:
```
# x y z radius
x1 y1 z1 r
x2 y2 z2 r
...
```

### Visualization

Convert the output to VTP format for visualization in ParaView:

```bash
cd python
python convert_output.py ../build/output.txt
```

The script supports the following options:
- `input_file`: (required) Input file containing particle positions and radii
- `--output-vtp`: Output VTP file name (default: `particles.vtp`)
- `--output-lethe`: Output Lethe insertion file name (default: `insertion_file.dat`)

Example with custom output names:
```bash
python convert_output.py output.txt --output-vtp my_particles.vtp --output-lethe my_insertion.dat
```

This will generate:
- `particles.vtp` (or custom name): VTK PolyData file for ParaView visualization
- `insertion_file.dat` (or custom name): Lethe-compatible particle insertion file

## References

1. Bitzek, E., Koskinen, P., Gähler, F., Moseler, M., & Gumbsch, P. (2006). *Structural Relaxation Made Simple*. Physical Review Letters, 97(17), 170201. [DOI: 10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201)

2. Guénolé, J., Nöhring, W. G., Vaid, A., Houllé, F., Xie, Z., Prakash, A., & Bitzek, E. (2020). *Assessment and optimization of the fast inertial relaxation engine (fire) for energy minimization in atomistic simulations and its implementation in lammps*. Computational Materials Science, 175, 109584.

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.
