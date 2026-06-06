/**
 * @file BlockJacobi_solver.hpp
 * @brief Definition of the BlockJacobi_solver class for solving the Poisson equation using Block Jacobi.
 */
#ifndef BLOCKJACOBI_SOLVER_HPP
#define BLOCKJACOBI_SOLVER_HPP

#include "Jacobi_solver.hpp"
#include <vector>

/**
 * @class BlockJacobi_solver
 * @brief A class to solve the Poisson equation using the Block Jacobi method.
 * 
 * This class inherits from Jacobi_solver and overrides the solve method to 
 * perform a Block Jacobi iteration, utilizing Eigen for exact local block solves.
 */
class BlockJacobi_solver : public Jacobi_solver {
public:
    /**
     * @brief Constructor for the Block Jacobi solver.
     * 
     * @param tol Convergence tolerance.
     * @param max_iter Maximum number of iterations.
     * @param n Number of grid points along each dimension.
     * @param source The right-hand side source function f(x,y).
     * @param boundary The boundary condition function g(x,y), defaulting to 0.0.
     */
    BlockJacobi_solver(double tol, unsigned max_iter, unsigned n, std::function<double(double,double)> source, std::function<double(double,double)> boundary = [](double, double){ return 0.0; }) 
        : Jacobi_solver(tol, max_iter, n, source, boundary) {};

    /**
     * @brief Solves the problem iteratively using the Block Jacobi approach.
     * 
     * @return std::vector<double> The computed solution vector flattened in 1D.
     */
    std::vector<double> solve() override;
};

#endif // BLOCKJACOBI_SOLVER_HPP
