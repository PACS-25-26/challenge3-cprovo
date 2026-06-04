#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include <functional>
#include <vector>

class Jacobi_solver {
public:
    Jacobi_solver(double tol, unsigned max_iter, unsigned n, std::function<double(double,double)> source) 
        : tol(tol), max_iter(max_iter), n(n), h(1.0 / (n - 1)), f(source) {};

    std::vector<double> solve();

private:
    double tol;
    unsigned max_iter;
    unsigned n;
    double h;
    std::function<double(double,double)> f; 
};

#endif // JACOBI_SOLVER_HPP
