#!/bin/bash

# Load Spack modules if running on the Polimi cluster
if [ -f "/software/spack/share/spack/setup-env.sh" ]; then
    echo "Loading cluster environment modules..."
    source /software/spack/share/spack/setup-env.sh
    spack load gcc@15.2.0
    spack load intel-oneapi-tbb@2022.3.0
    spack load openmpi@5.0.8
fi

# Compile the code
cd ..
make clean
make

if [ ! -f "jacobi_solver" ]; then
    echo "Error: jacobi_solver failed to compile. Exiting."
    exit 1
fi

# Move to data folder for results
cd test
mkdir -p data

# Grid sizes to test
SIZES=(16 32 64 128 256)
PROCS=(1 2 4)

echo "# Scalability Test Results" > benchmark.md
echo "Testing on $(uname -n)" >> benchmark.md
echo "---" >> benchmark.md

for p in "${PROCS[@]}"; do
    echo "Running with $p MPI processes..."
    echo "## MPI Processes: $p" >> benchmark.md
    echo "| Grid Size (n) | Time (s) | L2 Error |" >> benchmark.md
    echo "| ------------- | -------- | -------- |" >> benchmark.md
    
    for n in "${SIZES[@]}"; do
        echo "  Grid size: $n"
        
        # We enforce OMP_NUM_THREADS=2 for hybrid testing, you can change this
        export OMP_NUM_THREADS=2 
        
        # Run the solver and capture output
        OUTPUT=$(mpirun -np $p ../jacobi_solver $n)
        
        # Extract Time and Error from output for BlockJacobi_Homo
        TIME=$(echo "$OUTPUT" | grep "\[BlockJacobi_Homo\]" | awk '{print $4}')
        ERR=$(echo "$OUTPUT" | grep "\[BlockJacobi_Homo\]" | awk '{print $8}')
        
        echo "| $n | $TIME | $ERR |" >> benchmark.md
        
        # Save VTK output with specific names
        for vtk_file in ../Jacobi_Homo.vtk ../BlockJacobi_Homo.vtk ../BlockJacobi_NonHomo.vtk; do
            if [ -f "$vtk_file" ]; then
                base_name=$(basename "$vtk_file" .vtk)
                mv "$vtk_file" "data/${base_name}_${p}procs_${n}n.vtk"
            fi
        done
    done
    echo "" >> benchmark.md
done

echo "Tests completed! Results saved in benchmark.md and VTK files in data/"
