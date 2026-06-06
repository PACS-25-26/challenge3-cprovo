# Results Analysis

Based on the scalability test results available in `RESULT.md`, we can deduce the following observations about the performance and numerical correctness of the implemented solvers:

## 1. Numerical Correctness and Convergence
As the grid size ($n$) increases, the $L_2$ error consistently decreases, regardless of the number of MPI processes used. This confirms that the numerical method is consistent and correctly implemented across the distributed processes. The parallel implementation effectively reconstructs the exact solution for both the homogeneous and non-homogeneous Poisson problems.

## 2. Execution Time vs. Grid Size
For a fixed number of MPI processes, the execution time increases as the grid size grows. This behavior is expected because a larger grid involves more unknowns, requiring more computations per iteration and potentially more iterations to converge to the specified tolerance.

## 3. Parallel Scalability and Communication Overhead
When increasing the number of MPI processes, the expected behavior is a reduction in execution time due to distributed computation. However, the results show that:
- For small grid sizes (e.g., $16, 32, 64, 128$), using 2 or 4 MPI processes is often slower than using a single process.
- For $n=256$, using 1 process takes ~1.08s, while using 4 processes takes ~1.18s, which is relatively comparable, but still lacks strong scaling. The anomaly with 2 processes taking ~6.69s could be attributed to load balancing issues, cache contention, or latency on the test node.

**Conclusion:** 
The primary conclusion is that for small to moderate grid sizes, the overhead of MPI communication (passing ghost rows and synchronizing via `MPI_Allreduce`) outweighs the computational speedup. The computation per process is too small compared to the latency of sending and receiving data across the network. To observe significant parallel speedup (strong scaling), the problem size $n$ needs to be sufficiently large (e.g., $n \ge 512$ or $1024$) such that the computation-to-communication ratio is much higher.
