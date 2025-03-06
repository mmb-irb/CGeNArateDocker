import numpy as np
import pandas as pd
import os
import argparse

def main(sequence_length):
    MEAN = np.zeros((sequence_length, sequence_length))
    STDEV = np.zeros((sequence_length, sequence_length))
    MIN = np.zeros((sequence_length, sequence_length))
    MAX = np.zeros((sequence_length, sequence_length))

    conti = 1
    for i in range(sequence_length - 1):
        patata = i + 1
        contj = 1 + i
        for j in range(patata, sequence_length):
            myfile = f"{conti}-{contj}.dat"
            contj += 1
            if os.path.exists(myfile):
                #interact = pd.read_table(myfile, skiprows=1, header=None, delim_whitespace=True)  # Skip the header row and use whitespace as delimiter
                interact = pd.read_table(myfile, skiprows=1, header=None, sep=r'\s+')  # Skip the header row and use whitespace as delimiter
                #print(f"First few rows of {myfile}:")
                #print(interact.head())

                # Ensure that we are accessing the correct column (second column in this case)
                try:
                    values = pd.to_numeric(interact.iloc[:, 1], errors='coerce')
                    MEAN[i, j] = values.mean()
                    MEAN[j, i] = values.mean()
                    STDEV[i, j] = values.std()
                    STDEV[j, i] = values.std()
                    MIN[i, j] = values.min()
                    MIN[j, i] = values.min()
                    MAX[i, j] = values.max()
                    MAX[j, i] = values.max()
                except Exception as e:
                    print(f"Error processing {myfile}: {e}")
        conti += 1

    np.savetxt("distanceMean.contactMapMEAN.dat", MEAN)
    np.savetxt("distanceMean.contactMapMIN.dat", MIN)
    np.savetxt("distanceMean.contactMapMAX.dat", MAX)
    np.savetxt("distanceMean.contactMapSTDEV.dat", STDEV)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process sequence length.")
    parser.add_argument("sequence_length", type=int, help="Length of the sequence")
    args = parser.parse_args()
    main(args.sequence_length)

