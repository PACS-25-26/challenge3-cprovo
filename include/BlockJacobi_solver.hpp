#ifndef BLOCKJACOBI_SOLVER_HPP
#define BLOCKJACOBI_SOLVER_HPP

#include "Jacobi_solver.hpp"
#include <vector>

class BlockJacobi_solver : public Jacobi_solver {
public:
    BlockJacobi_solver(double tol, unsigned max_iter, unsigned n, std::function<double(double,double)> source, std::function<double(double,double)> boundary = [](double, double){ return 0.0; }) 
        : Jacobi_solver(tol, max_iter, n, source, boundary) {};

    std::vector<double> solve() override;
};

#endif // BLOCKJACOBI_SOLVER_HPP
