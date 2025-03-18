import argparse

def filter_pdb_chains(input_file, output_file, chains_to_keep):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            # Check if the line is an ATOM or HETATM record and belongs to the desired chains
            if line.startswith(("ATOM", "HETATM")) and line[21] in chains_to_keep:
                outfile.write(line)
            # Retain header and footer information like REMARK, HEADER, END
            elif not line.startswith(("ATOM", "HETATM")):
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description="Filter PDB file to keep only specified chains.")
    parser.add_argument("input_file", type=str, help="Path to the input PDB file.")
    parser.add_argument("output_file", type=str, help="Path to the output PDB file.")
    parser.add_argument(
        "chains_to_keep",
        type=str,
        help="Chains to keep, specified as a comma-separated list (e.g., 'A,B')."
    )

    args = parser.parse_args()
    chains = set(args.chains_to_keep.split(","))

    filter_pdb_chains(args.input_file, args.output_file, chains)

if __name__ == "__main__":
    main()

