# Testing and Scalability

This directory contains scripts to evaluate the performance and scalability of the matrix-free parallel Jacobi solver.

## Prerequisites
- A working MPI installation (e.g., OpenMPI or MPICH).
- A C++17 compliant compiler with OpenMP support.
- **Python 3** with `matplotlib` for the plotting pipeline.

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

## How to Plot the Results
After generating the `RESULT.md` file, you can automatically process the benchmarks and generate performance plots using the provided Python pipeline.
1. Python Requirements
The plotting script is written in Python 3 and relies on the matplotlib library for visualization.
2. Execution and Commands
To parse the data and build the charts, execute the dedicated Python script from the root of this directory:

```bash
   python3 plot_results.py
   ```

## Examining the output
- VTK files can be opened using **ParaView**.
- You can compare the times inside `RESULT.md` to analyze the speedup obtained through MPI and OpenMP hybrid parallelism.
- The script generates a high-resolution image named scalability_plot.png in the local folder.The generated image contains two log-log subplots:
1. Execution Time vs Grid Size: Used to analyze the execution scaling and speedup of the MPI parallel processes as the grid dimension increases ($n$).
2. $L_2$ Error vs Grid Size: Used to verify that the mathematical convergence rate of the Jacobi method remains consistent across different core counts and partition configurations.