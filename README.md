# Challenge 3: Parallel Jacobi Solvers

This repository contains the solution for Challenge 3, which involves implementing and optimizing parallel Jacobi and Block Jacobi solvers for the Poisson equation ($- \Delta u = f$) using a hybrid MPI + OpenMP approach.

## Repository Structure

The project is organized into the following directories and files:

- **`include/`**: Header files containing class definitions.
  - `Jacobi_solver.hpp`: Base class for the standard matrix-free Jacobi solver.
  - `BlockJacobi_solver.hpp`: Derived class implementing the Block Jacobi solver using `Eigen`.
- **`src/`**: Source files with the implementation of the solvers.
  - `Jacobi_solver.cpp`: Implementation of the matrix-free Jacobi iteration and MPI communication.
  - `BlockJacobi_solver.cpp`: Implementation of the Block Jacobi method with local Eigen factorizations.
  - `main.cpp`: The main driver program testing the solvers with homogeneous and non-homogeneous cases, generating VTK output.
- **`test/`**: Scripts and files related to testing and performance evaluation.
  - `test_scalability.sh`: Bash script to build the code and run scalability tests across varying grid sizes and MPI processes.
  - `RESULT.md`: Output file containing the execution times and $L_2$ errors recorded during tests.
  - `RESULTS_ANALYSIS.md`: Discussion and conclusions drawn from the test results.
  - `Testing and Scalability.md`: Instructions on how to run tests.
  - `test_jacobi.py` / `plot_results.py`: Python scripts for testing and visualization.
- **`Makefile`**: Build script to compile the project.
- **`Challenge25-26-3.pdf`**: The official challenge description.

## How to Run the Code

### Prerequisites
- A working MPI implementation (e.g., OpenMPI).
- A C++17 compliant compiler with OpenMP support.
- `Eigen3` library for sparse matrix operations (required by Block Jacobi).

### Building the Project
You can build the main executable by running `make` in the root directory:
```bash
make
```
This will compile the sources and place the object files in `obj/` and generate an executable.

### Execution
You can run the program using `mpirun` or `mpiexec`. By default, the program uses a grid size of $64 \times 64$.
```bash
mpirun -np 4 ./main
```
To specify a custom grid size $n$, pass it as an argument:
```bash
mpirun -np 4 ./main 128
```

### Testing and Scalability Instructions
To run the automated scalability tests, follow the instructions located in `test/Testing and Scalability.md`:
```bash
cd test
chmod +x test_scalability.sh
./test_scalability.sh
```
This will generate outputs and update `RESULT.md`. The VTK files produced can be visualized using ParaView.
