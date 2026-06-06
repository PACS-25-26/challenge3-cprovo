import matplotlib.pyplot as plt
import re
from collections import defaultdict
import os

def parse_result_md(file_path):
    results = defaultdict(lambda: {"n": [], "time": [], "error": []})
    current_procs = None
    
    with open(file_path, 'r') as f:
        for line in f:
            # Match process count
            proc_match = re.match(r'## MPI Processes:\s*(\d+)', line)
            if proc_match:
                current_procs = int(proc_match.group(1))
                continue
                
            # Match table row
            if current_procs and '|' in line:
                parts = [p.strip() for p in line.split('|')]
                if len(parts) >= 4 and parts[1].isdigit():
                    n = int(parts[1])
                    time = float(parts[2])
                    err = float(parts[3])
                    results[current_procs]["n"].append(n)
                    results[current_procs]["time"].append(time)
                    results[current_procs]["error"].append(err)
                    
    return results

def plot_scalability(results):
    if not results:
        print("No data found to plot. Run test_scalability.sh first.")
        return

    plt.figure(figsize=(12, 5))
    
    # Time vs N
    plt.subplot(1, 2, 1)
    for procs, data in sorted(results.items()):
        plt.plot(data["n"], data["time"], marker='o', label=f"{procs} Procs")
    plt.xlabel("Grid Size (n)")
    plt.ylabel("Time (s)")
    plt.title("Execution Time vs Grid Size")
    plt.legend()
    plt.grid(True)
    plt.xscale('log', base=2)
    plt.yscale('log')

    # L2 Error vs N
    plt.subplot(1, 2, 2)
    for procs, data in sorted(results.items()):
        plt.plot(data["n"], data["error"], marker='s', linestyle='--', label=f"{procs} Procs")
    plt.xlabel("Grid Size (n)")
    plt.ylabel("L2 Error")
    plt.title("L2 Error vs Grid Size")
    plt.legend()
    plt.grid(True)
    plt.xscale('log', base=2)
    plt.yscale('log')
    
    plt.tight_layout()
    plot_file = "scalability_plot.png"
    plt.savefig(plot_file)
    print(f"Plot saved successfully to {plot_file}")

if __name__ == "__main__":
    md_file = "RESULT.md"
    if os.path.exists(md_file):
        data = parse_result_md(md_file)
        plot_scalability(data)
    else:
        print(f"File {md_file} not found.")
