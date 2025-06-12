// Serial implementation of the 2D wave equation using finite difference (with timing)
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <chrono>

const double c = 1.0;               // wave speed
const double domain_size = 1.0;     // domain size [0,1] x [0,1]
const int Nx = 100;                 // number of grid points in x
const int Ny = 100;                 // number of grid points in y
const double dx = domain_size / (Nx - 1);
const double dy = domain_size / (Ny - 1);
const double dt = 0.4 * std::min(dx, dy) / c; // stability (Courant)
const int Nt = 500;                 // number of time steps
const double pi = 3.141592653589793;

double initial_condition(double x, double y) {
    double cx = 0.5, cy = 0.5, r = 0.1;
    return std::exp(-((x - cx)*(x - cx) + (y - cy)*(y - cy)) / (r * r));
}

inline int idx(int i, int j) {
    return i * Ny + j;
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<double> u_prev(Nx * Ny, 0.0);
    std::vector<double> u_curr(Nx * Ny, 0.0);
    std::vector<double> u_next(Nx * Ny, 0.0);

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double x = i * dx;
            double y = j * dy;
            u_prev[idx(i,j)] = initial_condition(x, y);
            u_curr[idx(i,j)] = u_prev[idx(i,j)];
        }
    }

    double cx2 = (c * dt / dx) * (c * dt / dx);
    double cy2 = (c * dt / dy) * (c * dt / dy);

    for (int n = 0; n < Nt; ++n) {
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                u_next[idx(i,j)] = 2.0 * u_curr[idx(i,j)] - u_prev[idx(i,j)] +
                    cx2 * (u_curr[idx(i+1,j)] - 2.0 * u_curr[idx(i,j)] + u_curr[idx(i-1,j)]) +
                    cy2 * (u_curr[idx(i,j+1)] - 2.0 * u_curr[idx(i,j)] + u_curr[idx(i,j-1)]);
            }
        }

        for (int i = 0; i < Nx; ++i) {
            u_next[idx(i,0)] = u_next[idx(i,Ny-1)] = 0.0;
        }
        for (int j = 0; j < Ny; ++j) {
            u_next[idx(0,j)] = u_next[idx(Nx-1,j)] = 0.0;
        }

        std::swap(u_prev, u_curr);
        std::swap(u_curr, u_next);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::ofstream fout("wave_output.dat");
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            fout << i * dx << " " << j * dy << " " << u_curr[idx(i,j)] << "\n";
        }
        fout << "\n";
    }
    fout.close();

    std::ofstream perf("performance_serial.txt", std::ios::app);
    perf << "Time: " << elapsed.count() << " seconds\n";
    perf.close();

    std::cout << "Simulation complete. Elapsed time: " << elapsed.count() << " seconds.\n";
    return 0;
}
