import sys

def process_pdb_file(input_file, output_file):
    dna_chain_ids = ["A", "B"]  # Chain IDs for the two DNA strands
    nn_chain_ids = ["C", "D", "E", "F", "G"]  # Chain IDs for NN residues
    current_chain = dna_chain_ids[0]
    new_lines = []
    last_residue_number = None
    dna_residues = []
    nn_chain_index = 0
    dna_part = True
    
    # Step 1: Collect DNA residue numbers
    with open(input_file, "r") as pdb:
        for line in pdb:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                residue_name = line[17:20].strip()
                residue_number = int(line[22:26])
                if residue_name != "NN":
                    dna_residues.append(residue_number)  # Collect DNA residue numbers
            
    # Calculate the midpoint for the DNA strand split
    if dna_residues:
        dna_midpoint = dna_residues[len(dna_residues) // 2]
    
    # Step 2: Process the file
    with open(input_file, "r") as pdb:
        for line in pdb:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                residue_name = line[17:20].strip()
                residue_number = int(line[22:26])
                
                if residue_name != "NN":  # DNA part
                    if residue_number == dna_midpoint and dna_part:
                        # Add TER tag and switch to the next DNA chain
                        new_lines.append(f"TER\n")
                        current_chain = dna_chain_ids[1]
                        dna_part = False
                    line = line[:21] + current_chain + line[22:]
                else:  # NN residue
                    if residue_number != last_residue_number:  # Add TER for new NN segment
                        new_lines.append(f"TER\n")
                        current_chain = nn_chain_ids[nn_chain_index]
                        nn_chain_index += 1
                    line = line[:21] + current_chain + line[22:]
                
                last_residue_number = residue_number
                new_lines.append(line)
            elif line.startswith("TER"):
                continue  # Skip existing TER lines
            else:
                new_lines.append(line)  # Preserve other lines (e.g., CRYST1)
    
    # Add a final TER tag at the end
    new_lines.append("TER\n")
    
    # Write the output file
    with open(output_file, "w") as out:
        out.writelines(new_lines)

# Main function to take input and output file parameters
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        process_pdb_file(input_file, output_file)
        print(f"Modified PDB file saved as {output_file}")

