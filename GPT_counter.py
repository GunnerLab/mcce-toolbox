import sys

def parse_pdb(file_path):
    """
    Parse a PDB file and count the number of amino acids, ligands, and waters.

    Args:
        file_path (str): Path to the PDB file.

    Returns:
        dict: Counts of amino acids, ligands, and waters.
    """
    amino_acids = set([
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    ])  # Standard amino acid residue names
    water_residue = 'HOH'  # Standard water residue name in PDB files

    amino_acid_count = 0
    ligand_count = 0
    water_count = 0

    try:
        with open(file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    # Extract the residue name (columns 18-20, 1-indexed)
                    residue_name = line[17:20].strip()

                    if line.startswith('ATOM') and residue_name in amino_acids:
                        # Count standard amino acids based on residue name
                        amino_acid_count += 1
                    elif line.startswith('HETATM'):
                        if residue_name == water_residue:
                            # Count water molecules explicitly
                            water_count += 1
                        else:
                            # Count other non-water HETATM entries as ligands
                            ligand_count += 1

        # Amino acids are counted based on standard residue names in ATOM records
        # Ligands are non-water, non-amino-acid HETATM entries
        # Waters are identified by the residue name 'HOH' in HETATM entries

        return {
            'amino_acids': amino_acid_count,
            'ligands': ligand_count,
            'waters': water_count
        }

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return None

# Example usage:
if __name__ == "__main__":
    try:
        input_path = sys.argv[1]  # PDB file path in terminal command
    except IndexError:
        print("\nPlease include the PDB file after the executable, e.g. 'ProtInfo_J.py 4lzt.pdb'.\n")
        raise
    counts = parse_pdb(input_path)
    if counts:
        print(f"Amino acids: {counts['amino_acids']}")
        print(f"Ligands: {counts['ligands']}")
        print(f"Waters: {counts['waters']}")

