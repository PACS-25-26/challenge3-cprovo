/**
 * @file Jacobi_solver.hpp
 * @brief Definition of the Jacobi_solver class for solving the Poisson equation.
 */
#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include <functional>
#include <vector>

/**
 * @class Jacobi_solver
 * @brief A class to solve the Poisson equation using the matrix-free parallel Jacobi method.
 * 
 * This base class implements a finite difference solver for the Poisson equation
 * over a 2D domain using a standard iterative Jacobi scheme with MPI and OpenMP parallelization.
 */
class Jacobi_solver {
public:
    /**
     * @brief Constructor for the Jacobi solver.
     * 
     * @param tol Convergence tolerance.
     * @param max_iter Maximum number of iterations.
     * @param n Number of grid points along each dimension.
     * @param source The right-hand side source function f(x,y).
     * @param boundary The boundary condition function g(x,y), defaulting to 0.0.
     */
    Jacobi_solver(double tol, unsigned max_iter, unsigned n, std::function<double(double,double)> source, std::function<double(double,double)> boundary = [](double, double){ return 0.0; }) 
        : tol(tol), max_iter(max_iter), n(n), h(1.0 / (n - 1)), f(source), g(boundary) {};

    /**
     * @brief Solves the problem iteratively.
     * 
     * @return std::vector<double> The computed solution vector flattened in 1D.
     */
    virtual std::vector<double> solve();
    
    /**
     * @brief Default virtual destructor.
     */
    virtual ~Jacobi_solver() = default;

protected:
    double tol;       ///< Tolerance for the stopping criterion.
    unsigned max_iter;///< Maximum allowed iterations.
    unsigned n;       ///< Grid size (number of points per dimension).
    double h;         ///< Grid spacing.
    std::function<double(double,double)> f; ///< Source function f(x, y).
    std::function<double(double,double)> g; ///< Boundary condition function g(x, y).
};

#endif // JACOBI_SOLVER_HPP
