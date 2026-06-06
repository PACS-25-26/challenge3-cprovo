#!/bin/bash

# Load Spack modules if running on the Polimi cluster
if [ -f "/software/spack/share/spack/setup-env.sh" ]; then
    echo "Loading cluster environment modules..."
    source /software/spack/share/spack/setup-env.sh
    spack load gcc@15.2.0
    spack load intel-oneapi-tbb@2022.3.0
    spack load openmpi@5.0.8
fi

# ==========================================
# GESTIONE PERCORSI SICURA
# ==========================================
# Trova la cartella reale in cui risiede questo script (cioè la cartella 'test')
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Trova la cartella principale del progetto (la root dove c'è il Makefile)
ROOT_DIR="$( cd "$SCRIPT_DIR/.." &> /dev/null && pwd )"

# Spostati nella root ed esegui il Make
cd "$ROOT_DIR"
echo "Compilazione in corso nella cartella: $(pwd)"

make clean
make

if [ ! -f "jacobi_solver" ]; then
    echo "Error: jacobi_solver failed to compile. Exiting."
    exit 1
fi

# Spostati nella cartella test per i risultati
cd "$SCRIPT_DIR"
mkdir -p data
# ==========================================

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
        
        export OMP_NUM_THREADS=2 
        
        # Esegui l'eseguibile puntando al percorso assoluto della root
        OUTPUT=$(mpirun -np $p "$ROOT_DIR/jacobi_solver" $n)
        
        TIME=$(echo "$OUTPUT" | grep "\[BlockJacobi_Homo\]" | awk '{print $4}')
        ERR=$(echo "$OUTPUT" | grep "\[BlockJacobi_Homo\]" | awk '{print $8}')
        
        echo "| $n | $TIME | $ERR |" >> benchmark.md
        
        # Salva i file VTK cercandoli nella cartella root dove vengono generati
        for vtk_file in "$ROOT_DIR"/Jacobi_Homo.vtk "$ROOT_DIR"/BlockJacobi_Homo.vtk "$ROOT_DIR"/BlockJacobi_NonHomo.vtk; do
            if [ -f "$vtk_file" ]; then
                base_name=$(basename "$vtk_file" .vtk)
                mv "$vtk_file" "data/${base_name}_${p}procs_${n}n.vtk"
            fi
        done
    done
    echo "" >> benchmark.md
done

echo "Tests completed! Results saved in benchmark.md and VTK files in data/"