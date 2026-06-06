/**
 * @file BlockJacobi_solver.cpp
 * @brief Implementation of the BlockJacobi_solver class.
 */
#include "BlockJacobi_solver.hpp"
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

/**
 * @brief Solves the Poisson equation using the Block Jacobi method.
 * 
 * In this implementation, the local domain block for each MPI process is solved
 * exactly using Eigen's SparseLU direct solver, while the global problem is 
 * iteratively solved using the Jacobi scheme for boundary updates between blocks.
 * 
 * @return std::vector<double> The global computed solution vector flattened in 1D.
 */
std::vector<double> BlockJacobi_solver::solve() {
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);   

    unsigned remainder = n % mpi_size;
    unsigned local_n = n / mpi_size;
    if ((unsigned)mpi_rank < remainder) {
        local_n += 1;
    }

    unsigned start_row = (n / mpi_size) * mpi_rank + std::min((unsigned)mpi_rank, remainder);

    std::vector<double> u((local_n + 2) * n, 0.0);

    // Initial guess and boundaries
    for (unsigned i = 1; i <= local_n; ++i) {
        unsigned global_i = start_row + i - 1;
        double x = global_i * h;
        for (unsigned j = 0; j < n; ++j) {
            double y = j * h;
            if (global_i == 0 || global_i == n - 1 || j == 0 || j == n - 1) {
                u[i * n + j] = g(x, y);
            }
        }
    }

    std::vector<double> f_val((local_n + 2) * n, 0.0);
    for (unsigned i = 1; i <= local_n; ++i) {
        unsigned global_i = start_row + i - 1;
        double x = global_i * h;
        for (unsigned j = 0; j < n; ++j) {
            double y = j * h;
            f_val[i * n + j] = f(x, y);
        }
    }

    // Build the local matrix A
    unsigned num_unknowns = local_n * n;
    Eigen::SparseMatrix<double> A(num_unknowns, num_unknowns);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(5 * num_unknowns);

    for (unsigned i = 1; i <= local_n; ++i) {
        unsigned global_i = start_row + i - 1;
        for (unsigned j = 0; j < n; ++j) {
            unsigned idx = (i - 1) * n + j;

            if (global_i == 0 || global_i == n - 1 || j == 0 || j == n - 1) {
                // Global boundary: Dirichlet condition (u = g)
                triplets.push_back({idx, idx, 1.0});
            } else {
                // Internal node
                triplets.push_back({idx, idx, 4.0});
                
                // Left neighbor
                triplets.push_back({idx, idx - 1, -1.0});
                
                // Right neighbor
                triplets.push_back({idx, idx + 1, -1.0});
                
                // Top neighbor (if i > 1, it's inside this block)
                if (i > 1) {
                    triplets.push_back({idx, idx - n, -1.0});
                }
                
                // Bottom neighbor (if i < local_n, it's inside this block)
                if (i < local_n) {
                    triplets.push_back({idx, idx + n, -1.0});
                }
            }
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    if(solver.info() != Eigen::Success) {
        std::cerr << "Rank " << mpi_rank << ": Factorization failed!\n";
    }

    double global_error = tol + 1.0;
    unsigned iter = 0;

    int top_neighbor = (mpi_rank > 0) ? mpi_rank - 1 : MPI_PROC_NULL;
    int bottom_neighbor = (mpi_rank < mpi_size - 1) ? mpi_rank + 1 : MPI_PROC_NULL;

    Eigen::VectorXd rhs(num_unknowns);
    Eigen::VectorXd x_sol(num_unknowns);

    while (global_error >= tol && iter < max_iter) {
        MPI_Sendrecv(&u[1 * n], n, MPI_DOUBLE, top_neighbor, 0,
                     &u[0 * n], n, MPI_DOUBLE, top_neighbor, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&u[local_n * n], n, MPI_DOUBLE, bottom_neighbor, 1,
                     &u[(local_n + 1) * n], n, MPI_DOUBLE, bottom_neighbor, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Build RHS
        #pragma omp parallel for
        for (unsigned i = 1; i <= local_n; ++i) {
            unsigned global_i = start_row + i - 1;
            for (unsigned j = 0; j < n; ++j) {
                unsigned idx = (i - 1) * n + j;
                
                if (global_i == 0 || global_i == n - 1 || j == 0 || j == n - 1) {
                    double x = global_i * h;
                    double y = j * h;
                    rhs(idx) = g(x, y);
                } else {
                    rhs(idx) = h * h * f_val[i * n + j];
                    
                    if (i == 1) {
                        rhs(idx) += u[0 * n + j]; 
                    }
                    if (i == local_n) {
                        rhs(idx) += u[(local_n + 1) * n + j]; 
                    }
                }
            }
        }

        x_sol = solver.solve(rhs);

        double local_error_sq = 0.0;
        
        #pragma omp parallel for reduction(+:local_error_sq)
        for (unsigned i = 1; i <= local_n; ++i) {
            unsigned global_i = start_row + i - 1;
            for (unsigned j = 0; j < n; ++j) {
                unsigned idx = (i - 1) * n + j;
                double new_val = x_sol(idx);
                
                if (!(global_i == 0 || global_i == n - 1 || j == 0 || j == n - 1)) {
                    double diff = new_val - u[i * n + j];
                    local_error_sq += diff * diff;
                }
                
                u[i * n + j] = new_val;
            }
        }

        double global_error_sq = 0.0;
        MPI_Allreduce(&local_error_sq, &global_error_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        global_error = std::sqrt(h * global_error_sq);
        iter++;
    }
    
    if (mpi_rank == 0) {
        std::cout << "Block Jacobi converged in " << iter << " iterations." << std::endl;
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
