#!/usr/bin/env python
"""
Created on Mar 13 09:00:00 2025

@author: Gehan Ranepura
"""

import os
import collections

# Directory containing the text files
directory = "pdb_output_mc_hbonds/"

# Dictionary to count occurrences of each hydrogen bond pair
hbond_counts = collections.defaultdict(int)
total_bonds = 0

# Iterate through all .txt files in the directory (excluding those starting with "blocking_file")
for filename in os.listdir(directory):
    if filename.endswith(".txt") and not filename.startswith("blocking_file"):
        filepath = os.path.join(directory, filename)
        with open(filepath, "r") as file:
            for line in file:
                parts = line.strip().split()
                if len(parts) >= 2:  # Ensure we have at least donor and acceptor
                    donor = parts[0]
                    acceptor = parts[1]
                    hbond = (donor, acceptor)
                    hbond_counts[hbond] += 1
                    total_bonds += 1

# Calculate percentages and write results to a file
output_file = "hbond_percentages.txt"
with open(output_file, "w") as out:
    out.write("Donor\t        Acceptor    \tPercentage of ms\n")
    for (donor, acceptor), count in sorted(hbond_counts.items(), key=lambda x: x[1], reverse=True):
        percentage = (count / total_bonds) * 100
        out.write(f"{donor}\t{acceptor}\t{percentage:.3f}%\n")

print(f"Hydrogen bond percentages saved to {output_file}.")

