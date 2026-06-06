#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include <functional>
#include <vector>

class Jacobi_solver {
public:
    Jacobi_solver(double tol, unsigned max_iter, unsigned n, std::function<double(double,double)> source, std::function<double(double,double)> boundary = [](double, double){ return 0.0; }) 
        : tol(tol), max_iter(max_iter), n(n), h(1.0 / (n - 1)), f(source), g(boundary) {};

    virtual std::vector<double> solve();
    virtual ~Jacobi_solver() = default;

protected:
    double tol;
    unsigned max_iter;
    unsigned n;
    double h;
    std::function<double(double,double)> f; 
    std::function<double(double,double)> g; 
};

#endif // JACOBI_SOLVER_HPP
