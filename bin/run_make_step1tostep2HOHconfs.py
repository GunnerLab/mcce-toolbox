#!/usr/bin/env python
"""
Created on Feb 11 01:25:00 2025

@author: Gehan Ranepura
"""

import os
import shutil
import subprocess

input_pdb = "step1_out.pdb"
output_pdb = "HOH_step2_out.pdb"
num_conformers = 25  # Change this if needed

# Initialize list for storing water molecule info
water_molecules = []

# Read input PDB file and extract water molecules (oxygen atoms)
with open(input_pdb) as file:
    for line in file:
        if line.startswith("HETATM") and line[17:20].strip() == "HOH" and line[12:16].strip() == "O":
            resi     = line[17:20]  # Residue
            chain_id = line[21]     # Chain ID
            res_num  = line[22:26]  # Residue number
            water_molecules.append((resi, chain_id, res_num))

# Ensure unique water molecules
water_molecules = list(set(water_molecules))

# Create a temporary directory for intermediate files
os.makedirs("temp_HOH_confs", exist_ok=True)
source_file = "HOH_confs.pdb"

# Run the script for each water molecule
for resi, chain_id, res_num in water_molecules:
    dir_name = f"temp_HOH_confs/{resi}_{chain_id}{res_num}"
    if os.path.exists(dir_name):
       shutil.rmtree(dir_name)  # Remove existing directory
    os.makedirs(dir_name)  # Create new directory

    temp_pdb = f"{dir_name}/{resi}_{chain_id}{res_num}.pdb"

    # Extract specific water molecule to a temporary file
    with open(temp_pdb, "w") as temp_file:
        with open(input_pdb) as file:
            for line in file:
                if (
                    line.startswith("HETATM")
                    and line[17:20].strip() == resi
                    and line[12:16].strip() == "O"
                    and line[21] == chain_id
                    and line[22:26].strip() == res_num
                ):
                    temp_file.write(line)

    # Call the function in the imported script
    subprocess.run(["make_step1tostep2HOHconfs.py", temp_pdb, "-N", str(num_conformers)])

    # Move the generated conformers into the temporary directory
    if os.path.exists(source_file):
        shutil.move(source_file, os.path.join(dir_name, "HOH_confs.pdb"))

# Merge all generated PDB files into a final output
with open(output_pdb, "w") as outfile:
    for resi, chain_id, res_num in water_molecules:
        dir_name = f"temp_HOH_confs/{resi}_{chain_id}{res_num}"
        confs_pdb_file = os.path.join(dir_name, "HOH_confs.pdb")
        if os.path.exists(confs_pdb_file):
            with open(confs_pdb_file) as infile:
                outfile.write(infile.read())

# Clean up the temporary files
shutil.rmtree("temp_HOH_confs")
print(f"Final output written to {output_pdb}")

