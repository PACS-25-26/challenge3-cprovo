/**
 * @file main.cpp
 * @brief Main entry point for the Poisson equation solvers evaluation.
 */
#include "Jacobi_solver.hpp"
#include "BlockJacobi_solver.hpp"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

const double PI = std::acos(-1.0);

// --- Homogeneous Problem ---

/**
 * @brief Forcing term f(x,y) for the homogeneous Dirichlet problem.
 * 
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @return double The value of the forcing term at (x,y).
 */
double f_source_homo(double x, double y) {
    return 8.0 * PI * PI * std::sin(2.0 * PI * x) * std::sin(2.0 * PI * y);
}

/**
 * @brief Exact solution u(x,y) for the homogeneous Dirichlet problem.
 * 
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @return double The value of the exact solution at (x,y).
 */
double u_exact_homo(double x, double y) {
    return std::sin(2.0 * PI * x) * std::sin(2.0 * PI * y);
}


// --- Non-Homogeneous Problem ---
// Equation: -Delta u = f
// Exact solution: u(x,y) = sin(2*pi*x)*sin(2*pi*y) + x + y
// Boundary conditions: g(x,y) = x + y
// Delta(x+y) = 0, so f is the same!

/**
 * @brief Forcing term f(x,y) for the non-homogeneous Dirichlet problem.
 * 
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @return double The value of the forcing term at (x,y).
 */
double f_source_nonhomo(double x, double y) {
    return 8.0 * PI * PI * std::sin(2.0 * PI * x) * std::sin(2.0 * PI * y);
}

/**
 * @brief Boundary condition function g(x,y) for the non-homogeneous Dirichlet problem.
 * 
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @return double The value of the boundary condition at (x,y).
 */
double g_boundary_nonhomo(double x, double y) {
    return x + y;
}

/**
 * @brief Exact solution u(x,y) for the non-homogeneous Dirichlet problem.
 * 
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @return double The value of the exact solution at (x,y).
 */
double u_exact_nonhomo(double x, double y) {
    return std::sin(2.0 * PI * x) * std::sin(2.0 * PI * y) + x + y;
}

/**
 * @brief Exports the solution in VTK format for visualization in ParaView.
 * 
 * @param filename Name of the output VTK file.
 * @param u The solution vector (1D flattened 2D grid).
 * @param n Grid size (number of points per dimension).
 * @param h Grid spacing.
 */
void export_vtk(const std::string& filename, const std::vector<double>& u, unsigned n, double h) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open file " << filename << " for writing.\n";
        return;
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "Laplace solution\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << n << " " << n << " 1\n";
    out << "ORIGIN 0.0 0.0 0.0\n";
    out << "SPACING " << h << " " << h << " 1.0\n";
    out << "POINT_DATA " << n * n << "\n";
    out << "SCALARS u double 1\n";
    out << "LOOKUP_TABLE default\n";

    for (unsigned j = 0; j < n; ++j) { // y-axis
        for (unsigned i = 0; i < n; ++i) { // x-axis
            out << u[i * n + j] << "\n";
        }
    }
    out.close();
}

/**
 * @brief Evaluates the computed solution against the exact solution.
 * 
 * Computes the L2 error and exports the solution to a VTK file (done by rank 0).
 * 
 * @param name The base name for output messages and files.
 * @param u The computed solution vector.
 * @param exact The exact solution function.
 * @param n Grid size.
 * @param mpi_rank The MPI rank of the current process.
 * @param time Elapsed time for the solver.
 */
void evaluate_solution(const std::string& name, const std::vector<double>& u, std::function<double(double,double)> exact, unsigned n, int mpi_rank, double time) {
    if (mpi_rank == 0) {
        double h = 1.0 / (n - 1);
        double l2_err_sq = 0.0;
        for (unsigned i = 0; i < n; ++i) {
            double x = i * h;
            for (unsigned j = 0; j < n; ++j) {
                double y = j * h;
                double diff = u[i * n + j] - exact(x, y);
                l2_err_sq += diff * diff;
            }
        }
        double l2_err = std::sqrt(h * l2_err_sq);
        std::cout << "[" << name << "] Time elapsed: " << time << " s, L2 Error: " << l2_err << std::endl;
        export_vtk(name + ".vtk", u, n, h);
    }
}

/**
 * @brief Main function.
 * 
 * Initializes MPI, parses command line arguments for grid size, instantiates solvers,
 * executes them, and triggers result evaluation.
 * 
 * @param argc Number of command-line arguments.
 * @param argv Command-line arguments.
 * @return int Exit status.
 */
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    unsigned n = 64; // Default size (reduced slightly to make Block Jacobi faster for testing)
    if (argc > 1) {
        n = std::atoi(argv[1]);
    }

    double tol = 1e-4;
    unsigned max_iter = 100000;

    if (mpi_rank == 0) {
        std::cout << "Grid size (n): " << n << "x" << n << std::endl;
        std::cout << "Tolerance: " << tol << ", Max Iterations: " << max_iter << std::endl;
        std::cout << "---" << std::endl;
    }

    // 1. Standard Jacobi - Homogeneous
    {
        Jacobi_solver solver(tol, max_iter, n, f_source_homo);
        double t_start = MPI_Wtime();
        std::vector<double> u = solver.solve();
        double t_end = MPI_Wtime();
        evaluate_solution("Jacobi_Homo", u, u_exact_homo, n, mpi_rank, t_end - t_start);
    }

    // 2. Block Jacobi - Homogeneous
    {
        BlockJacobi_solver solver(tol, max_iter, n, f_source_homo);
        double t_start = MPI_Wtime();
        std::vector<double> u = solver.solve();
        double t_end = MPI_Wtime();
        evaluate_solution("BlockJacobi_Homo", u, u_exact_homo, n, mpi_rank, t_end - t_start);
    }

    // 3. Block Jacobi - Non-Homogeneous
    {
        BlockJacobi_solver solver(tol, max_iter, n, f_source_nonhomo, g_boundary_nonhomo);
        double t_start = MPI_Wtime();
        std::vector<double> u = solver.solve();
        double t_end = MPI_Wtime();
        evaluate_solution("BlockJacobi_NonHomo", u, u_exact_nonhomo, n, mpi_rank, t_end - t_start);
    }

    MPI_Finalize();
    return 0;
}