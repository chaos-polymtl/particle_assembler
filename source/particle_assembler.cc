#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

// ----------------------------
// Parameters
// ----------------------------
const int N = 70486; // number of particles
const double phi_target = 0.30;
const double r_in = 0.0064;
const double thickness = 0.0238 - r_in;
const double height = 0.25;

const double k_pair = 1e3;
const double k_wall = 1e3;
const double grow_rate = 1.02;
const int fire_max_steps = 100000;
const double fire_dt = 1e-5;

// FIRE parameters (tunable)
const double FIRE_DT_INIT = fire_dt;
const double FIRE_DT_MAX = fire_dt * 200.0;
const double FIRE_ALPHA_INIT = 0.1;
const double FIRE_FINC = 1.1;
const double FIRE_FDEC = 0.5;
const double FIRE_ALPHA_DEC = 0.99;
const int FIRE_NMIN = 5;
const double FIRE_FTOL = 1e-6; // convergence on max force

// ----------------------------
// Particle structure
// ----------------------------
struct Particle {
  double x, y, z;    // position
  double vx, vy, vz; // velocity
  double fx, fy, fz; // force
};

// ----------------------------
// Utility functions
// ----------------------------
inline double sqr(double x) { return x * x; }

double cylinder_shell_volume(double r_in, double thickness, double height) {
  double r_out = r_in + thickness;
  return M_PI * (r_out * r_out - r_in * r_in) * height;
}

double radius_for_phi(int N, double phi_target, double r_in, double thickness,
                      double height) {
  double vol_shell = cylinder_shell_volume(r_in, thickness, height);
  double r3 = (phi_target * vol_shell) / (N * (4.0 / 3.0) * M_PI);
  return std::cbrt(r3);
}

// ----------------------------
// Linked cell grid
// ----------------------------
struct CellGrid {
  int nx, ny, nz;
  double cell_size;
  double x_min, y_min, z_min;
  std::vector<std::vector<int>> cells;

  CellGrid()
      : nx(0), ny(0), nz(0), cell_size(0.0), x_min(0.0), y_min(0.0),
        z_min(0.0) {}

  CellGrid(double r_out, double height, double cutoff) {
    init(r_out, height, cutoff);
  }

  void init(double r_out, double height, double cutoff) {
    // a conservative bounding box centered on origin in x,y and z=[0,height]
    cell_size = std::max(1e-12, cutoff); // avoid zero
    x_min = -r_out;
    y_min = -r_out;
    z_min = 0.0;
    nx = std::max(1, int(2 * r_out / cell_size) + 1);
    ny = std::max(1, int(2 * r_out / cell_size) + 1);
    nz = std::max(1, int(height / cell_size) + 1);
    cells.clear();
    cells.resize(nx * ny * nz);
  }

  inline int index(int ix, int iy, int iz) const {
    return (iz * ny + iy) * nx + ix;
  }

  void clear() {
    for (auto &c : cells)
      c.clear();
  }

  void insert(int p, const Particle &part) {
    int ix = int((part.x - x_min) / cell_size);
    int iy = int((part.y - y_min) / cell_size);
    int iz = int((part.z - z_min) / cell_size);
    if (ix >= 0 && ix < nx && iy >= 0 && iy < ny && iz >= 0 && iz < nz)
      cells[index(ix, iy, iz)].push_back(p);
    // if a particle falls outside bounds, we silently ignore insertion;
    // compute_forces will ignore pairs for particles outside cells, but
    // those particles will still have wall forces computed later.
  }
};

// ----------------------------
// Force calculation
// ----------------------------
// Computes pair forces (harmonic) and wall forces. Returns potential energy.
// Also computes and returns max_force via reference.
double compute_forces(std::vector<Particle> &parts, double r_particle,
                      double r_in, double thickness, double height,
                      CellGrid &grid, double &out_max_force) {
  int NN = (int)parts.size();
  for (auto &p : parts) {
    p.fx = p.fy = p.fz = 0.0;
  }

  grid.clear();
  for (int i = 0; i < NN; i++)
    grid.insert(i, parts[i]);

  double energy = 0.0;
  double r_out = r_in + thickness;
  double cutoff = 2.0 * r_particle * 1.1; // small skin

  // neighbor loops
  for (int ix = 0; ix < grid.nx; ix++) {
    for (int iy = 0; iy < grid.ny; iy++) {
      for (int iz = 0; iz < grid.nz; iz++) {
        int cell_index = grid.index(ix, iy, iz);
        auto &cell = grid.cells[cell_index];
        if (cell.empty())
          continue;

        // neighbor cells
        for (int dx = -1; dx <= 1; dx++) {
          for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
              int jx = ix + dx, jy = iy + dy, jz = iz + dz;
              if (jx < 0 || jy < 0 || jz < 0 || jx >= grid.nx ||
                  jy >= grid.ny || jz >= grid.nz)
                continue;
              int neigh_index = grid.index(jx, jy, jz);
              auto &neigh = grid.cells[neigh_index];

              for (int a : cell) {
                for (int b : neigh) {
                  if (b <= a)
                    continue; // avoid double counting
                  double dx = parts[a].x - parts[b].x;
                  double dy = parts[a].y - parts[b].y;
                  double dz = parts[a].z - parts[b].z;
                  double d2 = dx * dx + dy * dy + dz * dz;
                  if (d2 < cutoff * cutoff) {
                    double d = std::sqrt(d2);
                    // avoid division by zero
                    if (d < 1e-12) {
                      // tiny random jitter to separate coincident particles
                      static std::mt19937 rng(12345);
                      static std::uniform_real_distribution<double> udist(-1e-8,
                                                                          1e-8);
                      dx += udist(rng);
                      dy += udist(rng);
                      dz += udist(rng);
                      d = std::sqrt(dx * dx + dy * dy + dz * dz);
                      if (d < 1e-12)
                        d = 1e-12;
                    }
                    double overlap = 2.0 * r_particle - d;
                    if (overlap > 0) {
                      double fmag = k_pair * overlap;
                      double fx = fmag * dx / d;
                      double fy = fmag * dy / d;
                      double fz = fmag * dz / d;
                      parts[a].fx += fx;
                      parts[a].fy += fy;
                      parts[a].fz += fz;
                      parts[b].fx -= fx;
                      parts[b].fy -= fy;
                      parts[b].fz -= fz;
                      energy += 0.5 * k_pair * overlap * overlap;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // walls: inner, outer, z caps
  double max_force = 0.0;
  for (int i = 0; i < NN; i++) {
    Particle &p = parts[i];
    double r = std::sqrt(p.x * p.x + p.y * p.y);

    // inner wall (if inside forbidden core)
    double overlap_inner = r_in - r;
    if (overlap_inner > 0) {
      double fmag = k_wall * overlap_inner;
      double invr = (r > 1e-12) ? 1.0 / r : 0.0;
      p.fx += fmag * (p.x * invr);
      p.fy += fmag * (p.y * invr);
      energy += 0.5 * k_wall * overlap_inner * overlap_inner;
    }
    // outer wall (if outside)
    double overlap_outer = r - r_out;
    if (overlap_outer > 0) {
      double fmag = k_wall * overlap_outer;
      double invr = (r > 1e-12) ? 1.0 / r : 0.0;
      p.fx -= fmag * (p.x * invr);
      p.fy -= fmag * (p.y * invr);
      energy += 0.5 * k_wall * overlap_outer * overlap_outer;
    }
    // z caps
    if (p.z < 0) {
      double o = -p.z;
      p.fz += k_wall * o;
      energy += 0.5 * k_wall * o * o;
    }
    if (p.z > height) {
      double o = p.z - height;
      p.fz -= k_wall * o;
      energy += 0.5 * k_wall * o * o;
    }

    double fmag = std::sqrt(p.fx * p.fx + p.fy * p.fy + p.fz * p.fz);
    if (fmag > max_force)
      max_force = fmag;
  }

  out_max_force = max_force;
  return energy;
}

// ----------------------------
// FIRE minimizer
// ----------------------------
bool fire_minimize(std::vector<Particle> &parts, double r_particle, double r_in,
                   double thickness, double height, CellGrid &grid,
                   int max_steps, double dt_init, double ftol,
                   bool verbose = true) {
  const int NN = (int)parts.size();

  double dt = dt_init;
  double dt_max = std::max(dt * 10.0, FIRE_DT_MAX);
  double alpha = FIRE_ALPHA_INIT;
  double alpha_start = FIRE_ALPHA_INIT;
  int n_positive = 0;

  // initial compute (fills forces)
  double maxF = 0.0;
  double energy =
      compute_forces(parts, r_particle, r_in, thickness, height, grid, maxF);

  if (verbose) {
    std::cout << "  FIRE init: energy=" << energy << " maxF=" << maxF << "\n";
  }

  for (int step = 0; step < max_steps; ++step) {
    // 1) velocity update: v += dt * f  (mass = 1)
    for (int i = 0; i < NN; ++i) {
      parts[i].vx += parts[i].fx * dt;
      parts[i].vy += parts[i].fy * dt;
      parts[i].vz += parts[i].fz * dt;
    }

    // 2) compute P = sum v · f and norms
    double P = 0.0;
    double vnorm2 = 0.0;
    double fnorm2 = 0.0;
    for (int i = 0; i < NN; ++i) {
      P += parts[i].vx * parts[i].fx + parts[i].vy * parts[i].fy +
           parts[i].vz * parts[i].fz;
      vnorm2 += parts[i].vx * parts[i].vx + parts[i].vy * parts[i].vy +
                parts[i].vz * parts[i].vz;
      fnorm2 += parts[i].fx * parts[i].fx + parts[i].fy * parts[i].fy +
                parts[i].fz * parts[i].fz;
    }
    double vnorm = std::sqrt(vnorm2);
    double fnorm = std::sqrt(fnorm2);

    // 3) adapt dt and alpha
    if (P > 0.0) {
      n_positive++;
      if (n_positive > FIRE_NMIN) {
        dt = std::min(dt * FIRE_FINC, dt_max);
        alpha *= FIRE_ALPHA_DEC;
      }
    } else {
      n_positive = 0;
      dt *= FIRE_FDEC;
      alpha = alpha_start;
      // zero velocities on negative power
      for (int i = 0; i < NN; ++i) {
        parts[i].vx = parts[i].vy = parts[i].vz = 0.0;
      }
    }

    // 4) mix velocities: v = (1-alpha) v + alpha * (f * (vnorm/fnorm))
    if (vnorm > 1e-16 && fnorm > 1e-16) {
      double fac = vnorm / fnorm;
      for (int i = 0; i < NN; ++i) {
        parts[i].vx = (1.0 - alpha) * parts[i].vx + alpha * parts[i].fx * fac;
        parts[i].vy = (1.0 - alpha) * parts[i].vy + alpha * parts[i].fy * fac;
        parts[i].vz = (1.0 - alpha) * parts[i].vz + alpha * parts[i].fz * fac;
      }
    }

    // 5) integrate positions: r += v * dt
    for (int i = 0; i < NN; ++i) {
      parts[i].x += parts[i].vx * dt;
      parts[i].y += parts[i].vy * dt;
      parts[i].z += parts[i].vz * dt;
    }

    // 6) recompute forces at new positions
    energy =
        compute_forces(parts, r_particle, r_in, thickness, height, grid, maxF);

    // 7) convergence check
    if ((step % 100) == 0 && verbose) {
      std::cout << "    FIRE step " << step << " energy=" << energy
                << " maxF=" << maxF << " dt=" << dt << " alpha=" << alpha
                << "\n";
    }
    if (maxF < ftol) {
      if (verbose) {
        std::cout << "    FIRE converged in " << step << " steps (maxF=" << maxF
                  << ")\n";
      }
      return true; // converged
    }
  }

  // not converged
  if (verbose) {
    std::cout << "    FIRE did NOT converge after " << max_steps
              << " steps (final maxF unknown)\n";
  }
  return false;
}

int write_output_file(const std::vector<Particle> &parts, const double radius,
                      const std::string &filename) {
  // Open file for writing
  std::ofstream fout("output.txt");
  if (!fout) {
    std::cerr << "Error opening file!" << std::endl;
    return 1;
  }

  // Write header (optional)
  fout << "# x y z radius\n";

  // Write data in 3-column format
  for (unsigned int i = 0; i < parts.size(); ++i)
    fout << parts[i].x << " " << parts[i].y << " " << parts[i].z << " "
         << radius << "\n";

  fout.close();
  return 0;
}

// ----------------------------
// Main
// ----------------------------
int main() {
  double r_init = radius_for_phi(N, 0.02, r_in, thickness, height);
  double r_target = radius_for_phi(N, phi_target, r_in, thickness, height);
  double r_particle = r_init;

  std::cout << "Initial radius = " << r_init << " target radius = " << r_target
            << std::endl;

  std::mt19937 gen(42);
  std::uniform_real_distribution<double> u01(0.0, 1.0);

  std::vector<Particle> parts(N);
  double r_out = r_in + thickness;

  // random init
  for (int i = 0; i < N; i++) {
    double u = u01(gen);
    double radial = std::sqrt(u * (r_out * r_out - r_in * r_in) + r_in * r_in);
    double theta = u01(gen) * 2 * M_PI;
    double z = u01(gen) * height;
    parts[i].x = radial * std::cos(theta);
    parts[i].y = radial * std::sin(theta);
    parts[i].z = z;
    parts[i].vx = parts[i].vy = parts[i].vz = 0.0;
    parts[i].fx = parts[i].fy = parts[i].fz = 0.0;
  }

  // Growth + relaxation loop: after each growth step run FIRE to remove
  // overlaps
  int cycle = 0;
  while (r_particle < r_target) {
    cycle++;
    r_particle = std::min(r_particle * grow_rate, r_target);

    // rebuild grid using current cutoff based on the (grown) particle size
    double cutoff = 2.0 * r_particle * 1.1;
    CellGrid grid;
    grid.init(r_out, height, cutoff);

    std::cout << "\n=== Growth cycle " << cycle << ": r -> " << r_particle
              << " (cutoff=" << cutoff << ") ===\n";

    // Run FIRE to relax overlaps at this radius
    bool ok = fire_minimize(parts, r_particle, r_in, thickness, height, grid,
                            fire_max_steps, FIRE_DT_INIT, FIRE_FTOL,
                            /*verbose=*/true);

    // report energy after relax (one last force eval to get energy)
    double maxF = 0.0;
    double energy =
        compute_forces(parts, r_particle, r_in, thickness, height, grid, maxF);
    std::cout << "[cycle " << cycle << "] phi = "
              << (N * (4.0 / 3.0) * M_PI * r_particle * r_particle *
                  r_particle) /
                     cylinder_shell_volume(r_in, thickness, height)
              << ", energy = " << energy << ", maxF = " << maxF << "\n";

    // if FIRE didn't converge, you can try smaller grow_rate or more steps
    if (!ok) {
      std::cout
          << "Warning: FIRE failed to converge at this growth step. Consider "
             "reducing grow_rate or increasing fire_max_steps.\n";
      // we continue anyway (optionally you could break here)
    }

    write_output_file(parts, r_particle, "output.dat");
  }

  std::cout << "Packing finished.\n";
  return 0;
}
