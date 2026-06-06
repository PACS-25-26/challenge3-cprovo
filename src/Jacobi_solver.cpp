#include "Jacobi_solver.hpp"
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <cmath>

std::vector<double> Jacobi_solver::solve() {
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);   

    // Distribute rows among processes
    unsigned remainder = n % mpi_size;
    unsigned local_n = n / mpi_size;
    if ((unsigned)mpi_rank < remainder) {
        local_n += 1;
    }

    // Global start row for the current process
    unsigned start_row = (n / mpi_size) * mpi_rank + std::min((unsigned)mpi_rank, remainder);

    // Local grid including 1 top ghost row and 1 bottom ghost row
    // Row 0 is the top ghost row
    // Rows 1 to local_n are the local rows
    // Row local_n + 1 is the bottom ghost row
    std::vector<double> u((local_n + 2) * n, 0.0);
    std::vector<double> u_new((local_n + 2) * n, 0.0);

    // Initial guess and boundaries are 0.0, so no need to set them up for Dirichlet homogeneous case.
    for (unsigned i = 1; i <= local_n; ++i) {
        unsigned global_i = start_row + i - 1;
        double x = global_i * h;
        for (unsigned j = 0; j < n; ++j) {
            double y = j * h;
            if (global_i == 0 || global_i == n - 1 || j == 0 || j == n - 1) {
                u[i * n + j] = g(x, y);
                u_new[i * n + j] = g(x, y);
            }
        }
    }

    // Precompute the forcing term f
    std::vector<double> f_val((local_n + 2) * n, 0.0);
    for (unsigned i = 1; i <= local_n; ++i) {
        unsigned global_i = start_row + i - 1;
        double x = global_i * h;
        for (unsigned j = 0; j < n; ++j) {
            double y = j * h;
            f_val[i * n + j] = f(x, y);
        }
    }

    double global_error = tol + 1.0;
    unsigned iter = 0;

    int top_neighbor = (mpi_rank > 0) ? mpi_rank - 1 : MPI_PROC_NULL;
    int bottom_neighbor = (mpi_rank < mpi_size - 1) ? mpi_rank + 1 : MPI_PROC_NULL;

    while (global_error >= tol && iter < max_iter) {
        // Exchange ghost rows
        // Send our top row (1) to top neighbor's bottom ghost row (local_n + 1)
        // Receive top neighbor's bottom row into our top ghost row (0)
        MPI_Sendrecv(&u[1 * n], n, MPI_DOUBLE, top_neighbor, 0,
                     &u[0 * n], n, MPI_DOUBLE, top_neighbor, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Send our bottom row (local_n) to bottom neighbor's top ghost row (0)
        // Receive bottom neighbor's top row into our bottom ghost row (local_n + 1)
        MPI_Sendrecv(&u[local_n * n], n, MPI_DOUBLE, bottom_neighbor, 1,
                     &u[(local_n + 1) * n], n, MPI_DOUBLE, bottom_neighbor, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        double local_error_sq = 0.0;

        #pragma omp parallel for reduction(+:local_error_sq)
        for (unsigned i = 1; i <= local_n; ++i) {
            unsigned global_i = start_row + i - 1;
            
            // Do not update global boundary rows (Dirichlet conditions)
            if (global_i == 0 || global_i == n - 1) {
                continue; 
            }

            for (unsigned j = 1; j < n - 1; ++j) {
                // Mathematically correct Jacobi update for -Delta u = f
                u_new[i * n + j] = 0.25 * (u[(i - 1) * n + j] + u[(i + 1) * n + j] + 
                                           u[i * n + j - 1] + u[i * n + j + 1] + 
                                           h * h * f_val[i * n + j]);
                
                double diff = u_new[i * n + j] - u[i * n + j];
                local_error_sq += diff * diff;
            }
        }

        double global_error_sq = 0.0;
        MPI_Allreduce(&local_error_sq, &global_error_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        global_error = std::sqrt(h * global_error_sq);

        // Update u for next iteration
        #pragma omp parallel for
        for (unsigned i = 1; i <= local_n; ++i) {
            unsigned global_i = start_row + i - 1;
            if (global_i > 0 && global_i < n - 1) {
                for (unsigned j = 1; j < n - 1; ++j) {
                    u[i * n + j] = u_new[i * n + j];
                }
            }
        }
        
        iter++;
    }

    // Gather results
    std::vector<int> counts(mpi_size);
    std::vector<int> displs(mpi_size);

    int recv_count = local_n * n;
    MPI_Gather(&recv_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < mpi_size; ++i) {
            displs[i] = displs[i - 1] + counts[i - 1];
        }
    }

    std::vector<double> global_u;
    if (mpi_rank == 0) {
        global_u.resize(n * n, 0.0);
    }

    MPI_Gatherv(&u[1 * n], local_n * n, MPI_DOUBLE,
                global_u.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    return global_u;
}
