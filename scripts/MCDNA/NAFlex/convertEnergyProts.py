import csv
import argparse
import statistics

def convert_file(input_file, output_file, stats_file, divisor):
    # Open the input file and read its contents
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Initialize the Elastic Energy and V columns
    elastic_energy = []
    elastic_energy_per_bp = []

    # Process each line
    for i, line in enumerate(lines):
        # Extract values from each line
        values = line.strip().split(',')

        # The second value is the "Elastic energy"
        ee = float(values[3]) * 0.00239006
        
        # Store Elastic energy and calculate "Elastic energy per bp"
        elastic_energy.append(ee)
        elastic_energy_per_bp.append(ee / divisor)

    # Write the output to the CSV file
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        
        # Write the header with V1, V2, V3, ..., without extra quotes
        header = [''] + [f'V{i+1}' for i in range(len(lines))]
        writer.writerow(header)
        
        # Write the "Elastic energy" row
        ee_row = ["Elastic energy"] + [round(ee, 6) for ee in elastic_energy]  # Using the entire elastic_energy list
        writer.writerow(ee_row)
        
        # Write the "Elastic energy per bp" row
        ee_per_bp_row = ["Elastic energy per bp"] + [round(ee_bp, 6) for ee_bp in elastic_energy_per_bp]  # Using the entire elastic_energy_per_bp list
        writer.writerow(ee_per_bp_row)

    # Calculate mean and standard deviation for Elastic energy and Elastic energy per bp
    mean_ee = statistics.mean(elastic_energy)
    stdev_ee = statistics.stdev(elastic_energy) if len(elastic_energy) > 1 else 0

    mean_ee_per_bp = statistics.mean(elastic_energy_per_bp)
    stdev_ee_per_bp = statistics.stdev(elastic_energy_per_bp) if len(elastic_energy_per_bp) > 1 else 0

    # Write the statistics to a new file
    with open(stats_file, 'w', newline='') as statsfile:
        writer = csv.writer(statsfile)
        writer.writerow(['', 'Mean', 'SD'])
        writer.writerow(["Elastic energy (kcal/mol)", round(mean_ee, 6), round(stdev_ee, 6)])
        writer.writerow(["Elastic energy per bp (kcal/mol)", round(mean_ee_per_bp, 6), round(stdev_ee_per_bp, 6)])

    print(f"Files '{output_file}' and '{stats_file}' have been created successfully.")

# Set up argparse to handle command-line arguments
def main():
    parser = argparse.ArgumentParser(description="Convert input file to desired output format and compute statistics.")
    
    # Arguments for input and output files
    parser.add_argument('input_file', type=str, help="Path to the input file")
    parser.add_argument('output_file', type=str, help="Path to the output file")
    parser.add_argument('stats_file', type=str, help="Path to the statistics file")
    parser.add_argument('divisor', type=float, help="Divisor for calculating 'Elastic energy per bp'")
    
    # Parse arguments
    args = parser.parse_args()

    # Run the conversion
    convert_file(args.input_file, args.output_file, args.stats_file, args.divisor)

if __name__ == "__main__":
    main()

