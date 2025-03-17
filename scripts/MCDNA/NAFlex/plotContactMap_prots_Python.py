import os
import sys
import numpy as np

if len(sys.argv) != 8:
    print("Usage: python script.py <title> <ini> <end> <ini2> <end2> <offset> <outFile>")
    sys.exit(0)

# Parse arguments
title, ini, end, ini2, end2, offset, out = sys.argv[1:]
ini, end, ini2, end2, offset = map(int, [ini, end, ini2, end2, offset])

print(f"Title: {title}")

# Calculate matrix dimensions
pos = (end - ini) // offset + 1
pos2 = (end2 - ini2) // offset + 1

# Initialize matrices
mean_matrix = np.full((pos, pos2), np.nan)
stdev_matrix = np.full((pos, pos2), np.nan)
min_matrix = np.full((pos, pos2), np.nan)
max_matrix = np.full((pos, pos2), np.nan)

# Process files and calculate statistics
for i in range(ini, end + 1):
    index_i = i - ini
    for j in range(ini2, end2 + 1):
        index_j = j - ini2
        filename = f"{i}-{j}.dat"
        #print(filename)
        if os.path.exists(filename):
            data = np.loadtxt(filename, usecols=[1])  # Assume column 2 contains the data
            mean_matrix[index_i, index_j] = np.mean(data)
            stdev_matrix[index_i, index_j] = np.std(data)
            min_matrix[index_i, index_j] = np.min(data)
            max_matrix[index_i, index_j] = np.max(data)

# Save matrices to files with single-line, space-separated format
np.savetxt(f"{out}MEAN.dat", mean_matrix, fmt="%.6f", delimiter=" ")
np.savetxt(f"{out}STDEV.dat", stdev_matrix, fmt="%.6f", delimiter=" ")
np.savetxt(f"{out}MIN.dat", min_matrix, fmt="%.6f", delimiter=" ")
np.savetxt(f"{out}MAX.dat", max_matrix, fmt="%.6f", delimiter=" ")

print("Data processing complete. Matrices have been saved in single-line format.")
