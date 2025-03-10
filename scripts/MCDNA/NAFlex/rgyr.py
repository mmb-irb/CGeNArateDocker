import os
import argparse
import pandas as pd

def main(top, input_traj, output_file):
    # Write instructions to the input file
    with open('ptraj.in', 'w') as file:
        file.write(f"trajin {input_traj}\n")
        file.write(f"radgyr :* out {output_file}.dat\n")
        file.write("go\n")
    
    # Construct the command
    command = f"cpptraj {top} < ptraj.in > ptraj.out"
    print(f"Executing command: {command}")
    
    # Execute the command
    os.system(command)

    # File path for the input data
    file_path = f"{output_file}.dat" 

    # Read the input data from the file
    with open(file_path, "r") as file:
        rows = file.readlines()

    # Extract column headers and data rows
    headers = rows[0].split()
    data = [row.split() for row in rows[1:]]

    output_data = {}  # Start with an empty dictionary

    # Add "V1" and other columns based on `data` contents
    output_data["V1"] = [int(row[0]) for row in data]
    for col_idx in range(1, len(headers) - 1):
        column = [float(row[col_idx]) for row in data]
        output_data[f"V{col_idx+1}"] = column

    print(output_data)

    # Prepare the header dynamically based on the keys in the output_data dictionary
    header = [""] + [f"V{i}" for i in range(1, len(output_data) + 1)]

    # Write to the CSV file
    output_file_path = f"{output_file}.csv"
    with open(output_file_path, "w") as csv_file:
        # Write the header
        csv_file.write(",".join(f'"{col}"' for col in header) + "\n")
    
        # Write the first row: Index and column indices
        first_row = ['"1"'] + [str(i) for i in range(1, len(output_data["V1"]) + 1)]
        csv_file.write(",".join(first_row) + "\n")
    
        # Write the second row: Index and corresponding values
        second_row = ['"2"'] + [str(value) for value in output_data["V2"]]
        csv_file.write(",".join(second_row) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run cpptraj with input and output file parameters.")
    parser.add_argument("--top", required=True, help="Path to the topology file")
    parser.add_argument("--input_traj", required=True, help="Name of the input trajectory")
    parser.add_argument("--output_file", required=True, help="Name of the output file to be generated")
    
    args = parser.parse_args()
    main(args.top, args.input_traj, args.output_file)

