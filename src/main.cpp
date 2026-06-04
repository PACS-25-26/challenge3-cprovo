#include "Jacobi_solver.hpp"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

const double PI = std::acos(-1.0);

// Forcing term f(x,y)
double f_source(double x, double y) {
    return 8.0 * PI * PI * std::sin(2.0 * PI * x) * std::sin(2.0 * PI * y);
}

// Exact solution u(x,y)
double u_exact(double x, double y) {
    return std::sin(2.0 * PI * x) * std::sin(2.0 * PI * y);
}

// Export the solution in VTK format for ParaView
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

    for (unsigned i = 0; i < n; ++i) { // y-axis
        for (unsigned j = 0; j < n; ++j) { // x-axis
            out << u[i * n + j] << "\n";
        }
    }
    out.close();
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    unsigned n = 128; // Default size
    if (argc > 1) {
        n = std::atoi(argv[1]);
    }

    // Set tolerance and max iterations
    double tol = 1e-4;
    unsigned max_iter = 100000;

    Jacobi_solver solver(tol, max_iter, n, f_source);
    
    double t_start = MPI_Wtime();
    std::vector<double> u = solver.solve();
    double t_end = MPI_Wtime();

    if (mpi_rank == 0) {
        std::cout << "Grid size (n): " << n << "x" << n << std::endl;
        std::cout << "Time elapsed: " << t_end - t_start << " s" << std::endl;

        // Compute L2 error against exact solution
        double h = 1.0 / (n - 1);
        double l2_err_sq = 0.0;
        for (unsigned i = 0; i < n; ++i) {
            double x = i * h;
            for (unsigned j = 0; j < n; ++j) {
                double y = j * h;
                double diff = u[i * n + j] - u_exact(x, y);
                l2_err_sq += diff * diff;
            }
        }
        double l2_err = std::sqrt(h * l2_err_sq);
        std::cout << "L2 Error w.r.t exact solution: " << l2_err << std::endl;

        export_vtk("solution.vtk", u, n, h);
        std::cout << "Solution exported to solution.vtk" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
