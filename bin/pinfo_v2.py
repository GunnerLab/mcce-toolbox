#!/usr/bin/env python

from collections import defaultdict
from pathlib import Path
import re
import sys
from typing import TextIO, Tuple, Union


def line_counter(file_data: TextIO, source_mcce: bool = False) -> Tuple:
    """Accepts a string of text and returns counts for # of waters and amino acids, counting them by the number of lines they appear in.

    Args:
        file_data (string): file data in text form

    Returns:
        A 3-tuple: Number of waters, number of amino acids, list of ligand names
    """
    file_data = file_data.splitlines()
    waters = 0
    aa = 0
    non_aa = 0
    ligand_names = []
    # Standard amino acid residue names:
    amino_acids = set([
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    ])

    for this_line in file_data:

        # residues and non-residues:
        if this_line.startswith("HETATM") or this_line.startswith("ATOM"):
            residue_name = this_line[17:20].strip()

            # FIX: Remove if not displayed or part of output
            if this_line.startswith("ATOM") and residue_name not in amino_acids:
                # Count non-standard amino acids based on residue name
                non_aa += 1
        
        # FIX: CA can be calcium
        if " CA " in this_line and not ("REMARK" in this_line or "ANISOU" in this_line):
            aa += 1

        # waters and ligands names:
        if not source_mcce: 
            if this_line.startswith("HETATM") and "HOH" not in this_line:
                ligand_names.append(this_line[17:29])
            if this_line.startswith("HETATM") and "HOH" in this_line:
                # or "ANISOU" in this_line): # do we count ANISOU? Can"t remember
                waters += 1
        else:
            if this_line.startswith("HETATM") and "HOH" not in this_line:
                ligand_names.append(this_line[17:20])
            if this_line.startswith("HETATM") and "HOH" in this_line and " O " in this_line:
                waters += 1

    # print("Non-amino acid residues are this many: " + str(non_aa))
    return waters, ligand_names, aa


def display_table(data: list, graph_name: str = "Total", use_borders: bool = False) -> None:
    """Takes in a 3 x 3 list of data and outputs a nicely formatted table and displays amino acid, ligand, and water counts for the total PDB file and its chains.

    Args:
        data (3 x 3 list): Path to the PDB file.
        graph_name (string): what the top-left of the graph says, change it for files with multiple chains
        use_borders (bool): whether we want a border drawn around the table, or not

    Returns:
        Nothing (prints a table).
    """
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
        separator = "+" + "+".join("-" * (width + 2) for width in col_widths) + "+"

        # Build the table
        table = [separator]
        for row in table_data:
            row_line = "|" + "|".join(f" {str(row[col]).ljust(col_widths[col])} " for col in range(len(row))) + "|"
            table.append(row_line)
            table.append(separator)

    else:
        # Build the table without borders
        table = []
        for row in table_data:
            row_line = "     ".join(f"{str(row[col]).ljust(col_widths[col])}" for col in range(len(row)))
            table.append(row_line)

    # Print the table
    print("\n".join(table))


def parse_pdb_chains(pdb_filepath: str) -> Union[Tuple[list, list], None]:
    """Identifies and extracts individual chains from a PDB file.

    Args:
        pdb_filepath (str): Path to the PDB file.

    Returns:
        A 2-tuple: the chains' line data, and a list of the chains' names or
        None if an error occurs.
    """
    chains = defaultdict(list)

    try:
        with open(pdb_filepath, "r") as file:
            for line in file:
                # Only process lines starting with "ATOM" or "HETATM"
                if line.startswith(("ATOM", "HETATM")):
                    chain_id = line[21]  # Chain identifier is in 0-indexed column 22
                    #if chain_id not in chains:
                    #    chains[chain_id] = []
                    chains[chain_id].append(line)

        # Convert the chains dictionary to a list of strings
        chain_data = ["".join(chains[chain_id]) for chain_id in chains]

        return chain_data, list(chains.keys())

    except FileNotFoundError:
        print(f"Error: File {pdb_filepath!s} not found.")
        return None

    except Exception as e:
        print(f"An error occurred: {e}")
        return None


# TODO:
# 1. Place following code into a "main" function;
# 2. Do checks on number of sys.argv under condition: if __name__ == "__main__": section;
try: 
    input_file = Path(sys.argv[1]) 
    if not input_file.suffix:
        # will work whether or not the user includes .pdb
        input_file = input_file + ".pdb"

except IndexError:
    sys.exit("\nPlease include the PDB file after the executable, e.g. 'protinfo_v2.py 4lzt.pdb'.\n")

# maybe include multiple model warning- MCCE expects 1 model pdb file

try:
    with open(input_file) as file:  # for how many waters, ligands, aa"s, were in protein file
        data = file.read()

    with open("step1_out.pdb") as file:  # for how many waters, ligands, aa"s, post-processing
        mcce_data = file.read()

    with open("run1.log") as file:  # helps identify where changes have occurred
        log_data = file.read()

except FileNotFoundError:
    print(f"\nProtInfo requires {input_file}, its associated head1.lst, run1.log, and step1_out.pdb in the current directory.\n")
    raise SystemExit()

print("\nOutputting protinfo to 'protinfo2_output.txt'.\n")

sys.stdout = open("pinfo2_output.txt", "wt")

# get HOH, ligand names, and amino acid counts
HOH_count, ligand_names, aa_count = line_counter(data) 
HOH_mcce_count, ligand_mcce_names, aa_mcce_count = line_counter(mcce_data, source_mcce=True)

print("For the input file " + input_file + ", we find the following:\n")

# should count non-standard amino acid residues

# format the data for the function to pick up (maybe make this part of an Object)
# FIX: The display_table function should be responsible for this conversion to string:
interior_data = [
    [str(aa_count), str(aa_mcce_count), str(abs(aa_count - aa_mcce_count))],
    [str(len(set(ligand_names))), str(len(set(ligand_mcce_names))), str(abs(len(set(ligand_names)) - len(set(ligand_mcce_names))))],
    [str(HOH_count), str(HOH_mcce_count), str(abs(HOH_count - HOH_mcce_count))],
]
display_table(interior_data)
chains, chain_labels = parse_pdb_chains(input_file)
print(f"\nThe PDB file contains {len(chains)} chain(s).")

if len(chains) > 1:
    mcce_chains, _ = parse_pdb_chains("step1_out.pdb")

    # FIX: Why set the enumeration start to 1 only to reset it to 0-based when needed?
    for i, chain in enumerate(chains, start=1):
        print("") # formatting assistance

        # FIX: Consider shorter names, eg: chain_waters, chain_ligands, chain_aas 
        HOH_chain_count, ligand_chain_names, aa_chain_count = line_counter(chain) 

        # FIX: Consider shorter names, eg: mc_chain_waters, mc_chain_ligands, mc_chain_aas 
        HOH_mcce_chain_count, ligand_chain_mcce_names, aa_mcce_chain_count = line_counter(mcce_chains[i - 1], source_mcce=True) # adjust to 0-indexed
        # water count not happening properly for different chains for mcce_chains?

        chain_data = [
            [str(aa_chain_count), str(aa_mcce_chain_count), str(abs(aa_chain_count - aa_mcce_chain_count))],
            [str(len(set(ligand_chain_names))), str(len(set(ligand_chain_mcce_names))), str(abs(len(set(ligand_chain_mcce_names)) - len(set(ligand_chain_names))))],
            [str(HOH_chain_count), str(HOH_mcce_chain_count), str(abs(HOH_chain_count - HOH_mcce_chain_count))],
        ]
        display_table(chain_data, graph_name="Chain " + str(chain_labels[i - 1]))  # again, adjust to 0-indexed

# appears to be 00always_needed.tpl, so give the pathway to it
print("\n'00always_needed.tpl' is an editable file that controls ligand stripping",
      "and can be found in the /MCCE4/param directory.",
      "\nCustom changes to cofactor stripping thresholds may be passed in the -u option",
      "of step1.py: 'step1.py -u H2O_SASCUTOFF=0.10'"
      )

# == REVIEW PAUSED ==

log_data = log_data.splitlines()
# FIX: Change to defaultdict(tuple):
rules = {}  # dictionary to contain rules and examples for how molecules are renamed
err_top_files = ""
how_many_atoms_change = 0
missing_atoms = 0
NTR_line = ""
CTR_line = ""

for line in log_data:
    # how to account for multiple TR changes
    if "Labeling" in line and "NTR" in line:
        NTR_line += line + "\n"
    if "Labeling" in line and "CTR" in line:
        CTR_line += line + "\n"

    match = re.match(r"\s*Renaming \"(.*?)\" to \"(.*?)\"", line)
    if match:
        how_many_atoms_change += 1
        original, renamed = match.groups()

        # Determine the rule by finding the difference between original and renamed
        rule = ""
        for i, (o_char, r_char) in enumerate(zip(original, renamed)):
            if o_char != r_char:
                rule += f"Position {i}: {o_char!s} -> {r_char!s}\n"

        if len(original) != len(renamed):
            rule += f"Length difference: Original({len(original)}) -> Renamed({len(renamed)})\n"

        # Normalize the rule description
        rule = rule.strip()
        # Store one example for each unique rule
        if rule not in rules:
            rules[rule] = (original, renamed)

    if "Error!" in line:
        err_top_files += line + "\n"

    # count missing atoms, but don"t count the line saying "Missing heavy atoms detected."
    if "Missing" in line and "detected" not in line:
        missing_atoms += 1

print("\nThese residues have been modified:")
print("\nTERMINI:\n")
print(NTR_line + CTR_line)
print(str(missing_atoms) + " missing atoms added. This number includes atoms relabeled as NTR, or CTR. See run1.log for full list.")

if not ligand_names: # if ligand name list is empty, skip ligands
    msg = ("\nNo ligands detected.",
           "\nA list of all atoms that are modified can be found in run1.log.",
           f"\n{str(how_many_atoms_change)} atoms changed."
           )
    sys.exit(msg)  # finish the program early

print("\nLIGANDS:",
      f"{input_file} LGDNAMES: {str(set(ligand_names))}",  # might actually be nice to have
      f"step1_out.pdb LGDNAMES: {str(set(ligand_mcce_names))}",
      sep="\n"
      )

if rules:
    print("      The rules for changes, and examples:\n")
    for rule, example in rules.items():
        # print(f"Rule:\n{rule}")
        print(f"\t  {example[0]} -> {example[1]}")
    print("\nA list of all atoms that are modified can be found in run1.log.\n")

print(str(how_many_atoms_change) + " atoms changed.")

# find shared ligands between list of ligands, ligands w/o topology files
ligand_sifter = set((" ".join(ligand_names)).split()) & set(err_top_files.split())
ligand_names = [lig[:3] for lig in ligand_names]
ligand_top = set(ligand_names) - ligand_sifter

if ligand_top:
    print("\nWe have topology files for these ligands:")
    print("TPLFOUND: " + str(ligand_top))

if err_top_files:
    print("\nWe do not have topology files for these ligands:")
    print("\nNOTPL: ", [line[-3:] for line in err_top_files.splitlines()], end="")
    # to be concise, only print the 3-char names of the residues

    # FIX: Missing link?
    print("\nYou can either: (1) remove them from the input pdb file; (2) run with",
          "all atoms with zero charge; (3) make a topology file using instructions in:")

print("\nThanks for using protinfo: It's free!")
