#!/bin/bash

# Compile the code
cd ..
make clean
make

# Move to data folder for results
cd test
mkdir -p data

# Grid sizes to test
SIZES=(16 32 64 128 256)
PROCS=(1 2 4)

echo "# Scalability Test Results" > RESULT.md
echo "Testing on $(uname -n)" >> RESULT.md
echo "---" >> RESULT.md

for p in "${PROCS[@]}"; do
    echo "Running with $p MPI processes..."
    echo "## MPI Processes: $p" >> RESULT.md
    echo "| Grid Size (n) | Time (s) | L2 Error |" >> RESULT.md
    echo "| ------------- | -------- | -------- |" >> RESULT.md
    
    for n in "${SIZES[@]}"; do
        echo "  Grid size: $n"
        
        # We enforce OMP_NUM_THREADS=2 for hybrid testing, you can change this
        export OMP_NUM_THREADS=2 
        
        # Run the solver and capture output
        OUTPUT=$(mpirun -np $p ../jacobi_solver $n)
        
        # Extract Time and Error from output for BlockJacobi_Homo
        TIME=$(echo "$OUTPUT" | grep "\[BlockJacobi_Homo\]" | awk '{print $5}')
        ERR=$(echo "$OUTPUT" | grep "\[BlockJacobi_Homo\]" | awk '{print $9}')
        
        echo "| $n | $TIME | $ERR |" >> RESULT.md
        
        # Save VTK output with specific names
        for vtk_file in ../Jacobi_Homo.vtk ../BlockJacobi_Homo.vtk ../BlockJacobi_NonHomo.vtk; do
            if [ -f "$vtk_file" ]; then
                base_name=$(basename "$vtk_file" .vtk)
                mv "$vtk_file" "data/${base_name}_${p}procs_${n}n.vtk"
            fi
        done
    done
    echo "" >> RESULT.md
done

echo "Tests completed! Results saved in RESULT.md and VTK files in data/"
