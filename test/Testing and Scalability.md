# Testing and Scalability

This directory contains scripts to evaluate the performance and scalability of the matrix-free parallel Jacobi solver.

## Prerequisites
- A working MPI installation (e.g., OpenMPI or MPICH).
- A C++17 compliant compiler with OpenMP support.

## How to run the tests
1. Ensure the script has execution permissions:
   ```bash
   chmod +x test_scalability.sh
   ```
2. Run the scalability script:
   ```bash
   ./test_scalability.sh
   ```

The script will compile the code, create a `data` folder to store the output VTK files for different configurations, and populate `RESULT.md` with the execution times and $L_2$ errors.

## Examining the output
- VTK files can be opened using **ParaView**.
- You can compare the times inside `RESULT.md` to analyze the speedup obtained through MPI and OpenMP hybrid parallelism.
