import csv
import argparse

# Define a function to subtract the factor from "Elastic energy per bp" values
def adjust_elastic_energy(file_path, output_path, factor):
    with open(file_path, 'r') as input_file:
        reader = csv.reader(input_file)
        rows = list(reader)

    # Process each row to find "Elastic energy per bp"
    for i, row in enumerate(rows):
        if row[0].strip() == "Elastic energy per bp":
            # Subtract the factor from each value in the line (excluding the first element)
            rows[i] = [row[0]] + [str(float(value) - factor) for value in row[1:]]
            break

    # Write back the modified rows to the output file
    with open(output_path, 'w', newline='') as output_file:
        writer = csv.writer(output_file)
        writer.writerows(rows)

# Parse command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adjust 'Elastic energy per bp' values in the input file.")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("output_file", help="Path to save the output file")
    parser.add_argument("factor", type=float, help="Factor to subtract from 'Elastic energy per bp' values")
    args = parser.parse_args()

    # Call the function with arguments
    adjust_elastic_energy(args.input_file, args.output_file, args.factor)
    print(f"File updated and saved as {args.output_file}.")

