// Parallel MPI implementation of the 2D wave equation using finite difference (with timing output)
#include <mpi.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

const double c = 1.0;
const double domain_size = 1.0;
const int Nx = 100;
const int Ny = 100;
const double dx = domain_size / (Nx - 1);
const double dy = domain_size / (Ny - 1);
const double dt = 0.4 * std::min(dx, dy) / c;
const int Nt = 500;

inline int idx(int i, int j, int ny) { return i * ny + j; }

double initial_condition(double x, double y) {
    double cx = 0.5, cy = 0.5, r = 0.1;
    return std::exp(-((x - cx)*(x - cx) + (y - cy)*(y - cy)) / (r * r));
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time = 0.0, end_time = 0.0;
    if (rank == 0) start_time = MPI_Wtime();

    int rows_per_proc = Nx / size;
    int remainder = Nx % size;

    int local_Nx = rows_per_proc + (rank < remainder ? 1 : 0);
    int offset = rank * rows_per_proc + std::min(rank, remainder);

    std::vector<double> u_prev((local_Nx + 2) * Ny, 0.0);
    std::vector<double> u_curr((local_Nx + 2) * Ny, 0.0);
    std::vector<double> u_next((local_Nx + 2) * Ny, 0.0);

    for (int i = 1; i <= local_Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double x = (offset + i - 1) * dx;
            double y = j * dy;
            u_prev[idx(i, j, Ny)] = initial_condition(x, y);
            u_curr[idx(i, j, Ny)] = u_prev[idx(i, j, Ny)];
        }
    }

    double cx2 = (c * dt / dx) * (c * dt / dx);
    double cy2 = (c * dt / dy) * (c * dt / dy);

    for (int n = 0; n < Nt; ++n) {
        std::vector<MPI_Request> reqs;
        MPI_Request req;

        if (rank > 0) {
            MPI_Isend(&u_curr[idx(1, 0, Ny)], Ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req);
            reqs.push_back(req);
            MPI_Irecv(&u_curr[idx(0, 0, Ny)], Ny, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &req);
            reqs.push_back(req);
        }
        if (rank < size - 1) {
            MPI_Isend(&u_curr[idx(local_Nx, 0, Ny)], Ny, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &req);
            reqs.push_back(req);
            MPI_Irecv(&u_curr[idx(local_Nx + 1, 0, Ny)], Ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &req);
            reqs.push_back(req);
        }

        if (!reqs.empty())
            MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

        for (int i = 1; i <= local_Nx; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                u_next[idx(i, j, Ny)] = 2.0 * u_curr[idx(i, j, Ny)] - u_prev[idx(i, j, Ny)] +
                    cx2 * (u_curr[idx(i + 1, j, Ny)] - 2.0 * u_curr[idx(i, j, Ny)] + u_curr[idx(i - 1, j, Ny)]) +
                    cy2 * (u_curr[idx(i, j + 1, Ny)] - 2.0 * u_curr[idx(i, j, Ny)] + u_curr[idx(i, j - 1, Ny)]);
            }
        }

        for (int i = 1; i <= local_Nx; ++i) {
            u_next[idx(i, 0, Ny)] = 0.0;
            u_next[idx(i, Ny - 1, Ny)] = 0.0;
        }
        if (rank == 0) {
            for (int j = 0; j < Ny; ++j) u_next[idx(1, j, Ny)] = 0.0;
        }
        if (rank == size - 1) {
            for (int j = 0; j < Ny; ++j) u_next[idx(local_Nx, j, Ny)] = 0.0;
        }

        std::swap(u_prev, u_curr);
        std::swap(u_curr, u_next);
    }

    std::ofstream fout("wave_rank" + std::to_string(rank) + ".dat");
    for (int i = 1; i <= local_Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            fout << (offset + i - 1) * dx << " " << j * dy << " " << u_curr[idx(i, j, Ny)] << "\n";
        }
        fout << "\n";
    }
    fout.close();

    if (rank == 0) {
        end_time = MPI_Wtime();
        double elapsed = end_time - start_time;
        std::ofstream perf("performance_parallel.txt", std::ios::app);
        perf << "Processes: " << size << ", Time: " << elapsed << " seconds\n";
        perf.close();
        std::cout << "Elapsed time: " << elapsed << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}
