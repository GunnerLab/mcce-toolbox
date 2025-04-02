#!/usr/bin/env python
"""
Created on Mar 13 09:00:00 2025

@author: Gehan Ranepura
"""

import os
import collections
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Compute hydrogen bond percentages from text files.")
parser.add_argument("dir", nargs="?", type=str, default="pdb_output_mc_hbonds", help="Directory containing hbond connection files (default: pdb_output_mc_hbonds)")
args = parser.parse_args()

# Validate directory existence
if not os.path.isdir(args.dir):
    print(f"Error: Directory '{args.dir}' not found.")
    exit(1)

# Dictionary to count occurrences of each hydrogen bond pair
hbond_counts = collections.Counter()
pdb_count = 0  # Number of PDBs processed (i.e., number of .txt files)

# Iterate through all .txt files in the directory (excluding those starting with "blocking_file")
for entry in os.scandir(args.dir):
    if entry.is_file() and entry.name.endswith(".txt") and not entry.name.endswith("blocking_file.txt"):
        pdb_count += 1
        with open(entry.path, "r") as file:
            unique_bonds = set()  # To track unique bonds in this PDB file
            for line in file:
                parts = line.strip().split()
                if len(parts) >= 2:  # Ensure we have at least donor and acceptor
                    donor, acceptor = parts[:2]
                    hbond = (donor, acceptor)
                    unique_bonds.add(hbond)  # Add unique bond for this PDB
            for bond in unique_bonds:
                hbond_counts[bond] += 1  # Count bond only once per PDB

# Handle case where no valid PDBs were found
if pdb_count == 0:
    print("No valid PDB files found. Exiting.")
    exit(1)
print(f"Number of PDBs = {pdb_count}")

# Calculate percentages and write results to a file
dir_name = os.path.basename(os.path.normpath(args.dir))
output_file = f"{dir_name}_percentages.txt"

with open(output_file, "w") as out:
    out.write(f"{'Donor':<20}{'Acceptor':<20}{'PDB Count':<15}{'Percentage of PDBs (%)'}\n")
    out.write("=" * 75 + "\n")
    for (donor, acceptor), count in sorted(hbond_counts.items(), key=lambda x: x[1], reverse=True):
        percentage = (count / pdb_count) * 100
        out.write(f"{donor:<20}{acceptor:<20}{count:<15}{percentage:.3f}%\n")

print(f"Hydrogen bond percentages saved to {output_file}")

