#!/usr/bin/env python
import sys
import re

def line_counter(file_data, source_mcce=False):  
    """Accepts a string of text and returns counts for # of waters and amino acids, counting them by the number of lines they appear in.

    Args:
        file_data (string): file data in text form

    Returns:
        list: A list of strings, each containing the data for one chain."""
    
    file_data_lines = file_data.splitlines()
    waters = 0
    aa = 0
    non_aa = 0
    ligand_names = []
    amino_acids = set([
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    ])  # Standard amino acid residue names

    for this_line in file_data_lines:

        if " CA " in this_line and this_line.startswith("ATOM"):
            aa += 1
        if source_mcce == False: 
            if this_line.startswith("HETATM") and "HOH" not in this_line:
                ligand_names.append(this_line[17:29])
            if "HOH" in this_line and "HETATM" in this_line:
                waters += 1
        else:
            if this_line.startswith("HETATM") and "HOH" not in this_line:
                ligand_names.append(this_line[17:20])
            if this_line.startswith('HETATM') and "HOH" in this_line and " O " in this_line:
                waters += 1
        if this_line.startswith('HETATM') or this_line.startswith('ATOM'):
            residue_name = this_line[17:20].strip()

            if this_line.startswith('ATOM') and residue_name not in amino_acids:
                # Count standard amino acids based on residue name
                non_aa += 1

    # print("Non-amino acid residues are this many: " + str(non_aa))

    return waters, ligand_names, aa

def display_table(data, graph_name="Total", use_borders=False):
    """Takes in a 3 x 3 list of data and outputs a nicely formatted table. Used to display amino acid, ligand, and water counts for the total PDB file and its chains.

    Args:
        data (3 x 3 list): Path to the PDB file.
        graph_name (string): what the top-left of the graph says, change it for files with multiple chains
        use_borders (bool): whether we want a border drawn around the table, or not

    Returns:
        list: A list of strings, each containing the data for one chain."""

    if len(data) != 3 or any(len(row) != 3 for row in data):
        raise ValueError("Input data must be a 3x3 matrix.")

    # Add headers and footers to create a 4x4 table
    table_data = [
        [graph_name, "PDB Count", "MCCE Count", "Difference"],  # Header row
        ["Amino Acids"] + data[0],  # First data row
        ["Ligands"] + data[1],  # Second data row
        ["Waters"] + data[2],  # Third data row
    ]

    # Determine column widths
    col_widths = [max(len(str(row[col])) for row in table_data) for col in range(len(table_data[0]))]
    if use_borders:
        # Create a separator line
        separator = '+' + '+'.join('-' * (width + 2) for width in col_widths) + '+'

        # Build the table
        table = [separator]
        for row in table_data:
            row_line = '|' + '|'.join(f" {str(row[col]).ljust(col_widths[col])} " for col in range(len(row))) + '|'
            table.append(row_line)
            table.append(separator)

    else:

        # Build the table without borders
        table = []
        for row in table_data:
            row_line = '     '.join(f"{str(row[col]).ljust(col_widths[col])}" for col in range(len(row)))
            table.append(row_line)

    # Print the table
    print('\n'.join(table))

def parse_pdb_chains(pdb_file):
    """Identifies and extracts individual chains from a PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        list: A list of strings, each containing the data for one chain."""
    
    chains = {}

    try:
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith(('ATOM', 'HETATM')):
                    chain_id = line[21]  # Chain identifier is in column 22 (0-indexed: 21)
                    if chain_id not in chains:
                        chains[chain_id] = []
                    chains[chain_id].append(line)

        # Convert the chains dictionary to a list of strings
        chain_data = ["".join(chains[chain_id]) for chain_id in chains]

        return chain_data, list(chains.keys()) # return the chains, and a list of the chains' names

    except FileNotFoundError:
        print(f"Error: File '{pdb_file}' not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []

if __name__ == "__main__":

    try:
        input_file = sys.argv[1].lower() # take the first terminal argument after the Python script is called. WHAT IF FILE NOT LOWERCASE? Include option to turn the .lower() off
        if ".pdb" not in input_file: # will work whether or not the user includes .pdb
            input_file = input_file + ".pdb"
    except IndexError:
        print("\nPlease include the PDB file after the executable, e.g. 'ProtInfo_J.py 4lzt.pdb'. protinfo expects lowercase file names.\n")
        raise SystemExit() # just quit the program

    # maybe include multiple model warning- MCCE expects 1 model pdb file

    try:
        with open(input_file, 'r') as file:  # for how many waters, ligands, aa's, were in protein file
            data = file.read()

        with open('step1_out.pdb', 'r') as file:  # for how many waters, ligands, aa's, post-processing
            mcce_data = file.read()

        with open('run1.log', 'r') as file:  # helps identify where changes have occurred
            log_data = file.read()

    except FileNotFoundError:
        print(f"\nProtInfo requires {input_file}, its associated, run1.log, and step1_out.pdb in the current directory.\n")
        raise SystemExit()

    # print("\nOutputting protinfo to 'protinfo2_output.txt'.\n")

    # get HOH, ligand names, and amino acid counts
    HOH_count, ligand_names, aa_count = line_counter(data) 
    HOH_mcce_count, ligand_mcce_names, aa_mcce_count = line_counter(mcce_data, source_mcce=True)

    print("\nNew ProtInfo by Jared! It's still in progress!\n")
    print("## For the input file " + input_file + ", we find the following:\n")

    # should count non-standard amino acid residues

    # format the data for the function to pick up (maybe make this part of an Object)
    interior_data = [
        [str(aa_count), str(aa_mcce_count), str(abs(aa_count - aa_mcce_count))],
        [str(len(set(ligand_names))), str(len(set(ligand_mcce_names))), str(abs(len(set(ligand_names)) - len(set(ligand_mcce_names))))],
        [str(HOH_count), str(HOH_mcce_count), str(abs(HOH_count - HOH_mcce_count))],
    ]

    display_table(interior_data)
    chains, chain_labels = parse_pdb_chains(input_file)
    print(f"\n## The PDB file contains {len(chains)} chain(s).")

    if len(chains) > 1:

        mcce_chains, _ = parse_pdb_chains('step1_out.pdb')

        for i, chain in enumerate(chains, start=1):
        
            print("") # formatting assistance

            HOH_chain_count, ligand_chain_names, aa_chain_count = line_counter(chain) 
            HOH_mcce_chain_count, ligand_chain_mcce_names, aa_mcce_chain_count = line_counter(mcce_chains[i - 1], source_mcce=True) # adjust to 0-indexed
            # water count not happening properly for different chains for mcce_chains?

            # there has to be a better way to do this
            chain_data = [
                [str(aa_chain_count), str(aa_mcce_chain_count), str(abs(aa_chain_count - aa_mcce_chain_count))],
                [str(len(set(ligand_chain_names))), str(len(set(ligand_chain_mcce_names))), str(abs(len(set(ligand_chain_mcce_names)) - len(set(ligand_chain_names))))],
                [str(HOH_chain_count), str(HOH_mcce_chain_count), str(abs(HOH_chain_count - HOH_mcce_chain_count))],
            ]
            display_table(chain_data, graph_name="Chain " + str(chain_labels[i - 1])) # again, adjust to 0-indexed

    print("\nThese residues are stripped if they are surface exposed. The percent exposure limit may be edited in '00always_needed.tpl'") # appears to be 00always_needed.tpl, so give the pathway to it

    log_data = log_data.splitlines()
    rules = {} # dictionary to contain rules and examples for how molecules are renamed
    err_top_files = ""
    how_many_atoms_change = 0
    missing_atoms = 0
    NTR_line = ""
    CTR_line = ""
    prox_ligands = ""

    for line in log_data:

        if "Distance below bond threshold" in line:
            prox_ligands += "   " + line + "\n"

        if "Labeling" in line and "NTR" in line:
            NTR_line += line + "\n" # how to account for multiple NTR changes
        if "Labeling" in line and "NTG" in line:
            NTR_line += line + " (NTR for GYL)" + "\n"
        if "Labeling" in line and "CTR" in line:
            CTR_line += line + "\n"

        match = re.match(r'\s*Renaming \"(.*?)\" to \"(.*?)\"', line)
        if match:
            how_many_atoms_change += 1
            original, renamed = match.groups()

            # Determine the rule by finding the difference between original and renamed
            rule = ""
            for i, (o_char, r_char) in enumerate(zip(original, renamed)):
                if o_char != r_char:
                    rule += f"Position {i}: '{o_char}' -> '{r_char}'\n"

            if len(original) != len(renamed):
                rule += f"Length difference: Original({len(original)}) -> Renamed({len(renamed)})\n"

            # Normalize the rule description
            rule = rule.strip()

            # Store one example for each unique rule
            if rule not in rules:
                rules[rule] = (original, renamed)

        if "Error!" in line:
            err_top_files += line + "\n" 

        # count missing atoms, but don't count the line saying "Missing heavy atoms detected."
        if "Missing" in line and "detected" not in line: 
            missing_atoms += 1

    print("\nThese residues have been modified:")
    print("\n## TERMINI:\n")
    print(NTR_line + CTR_line)
    print(str(missing_atoms) + " missing atoms added. This number includes atoms relabeled as NTR, or CTR. See run1.log for full list.")

    if not ligand_names: # if ligand name list is empty, skip ligands

        print("\nNo ligands detected.")
        print("\nA list of all atoms that are modified can be found in run1.log.")
        print("\n" + str(how_many_atoms_change) + " atoms changed.")
        raise SystemExit() # finish the program early

    print("\n## LIGANDS:")

    print(f"\n      {input_file} LGDNAMES: " + str(set(ligand_names)))
    print("\n      step1_out.pdb LGDNAMES: " + str(set(ligand_mcce_names)))
    print("\n" + prox_ligands)

    if bool(rules):

        print("      The rules for changes, and examples:\n")

        for rule, example in rules.items():
            # print(f"Rule:\n{rule}")
            print(f"      {example[0]} -> {example[1]}")

        print("\nA list of all atoms that are modified can be found in run1.log.\n")

    print(str(how_many_atoms_change) + " atoms changed.")

    # find shared ligands between list of ligands, ligands w/o topology files
    ligand_sifter = set((" ".join(ligand_names)).split()) & set(err_top_files.split()) 
    ligand_names = [l[0:3] for l in ligand_names] 
    ligand_top = set(ligand_names) - ligand_sifter

    if ligand_top:

        print("\nWe have topology files for these ligands:")
        print("      TPLFOUND: " + str(ligand_top))

    if err_top_files:

        print("\nWe do not have topology files for these ligands:")
        print("\nNOTPL: ", end='')
        print([line[-3:] for line in err_top_files.splitlines()]) # to be concise, only print the 3-char names of the residues
        print("\nYou can (1) remove them from the input pdb file; (2) run with all atoms with zero charge; (3) make a topology file using instructions in:")

    print("\nThanks for using protinfo!\n")
