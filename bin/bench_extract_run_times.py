#!/usr/bin/env python
"""
Created on Fri Dec 13 14:00:00 2024

@author: Gehan Ranepura
"""

import os
import pandas as pd
import argparse

def extract_run_times(runs_folder):
    # Initialize an empty list to store data
    data = []

    # Loop through each pdb folder in the runs directory
    for pdb_folder in os.listdir(runs_folder):
        pdb_folder_path = os.path.join(runs_folder, pdb_folder)

        # Check if it's a folder
        if os.path.isdir(pdb_folder_path):
            run_log_file = os.path.join(pdb_folder_path, 'run.log')
            rot_stat_file = os.path.join(pdb_folder_path, 'rot_stat')
            prot_pdb_file = os.path.join(pdb_folder_path, 'prot.pdb')

            # Skip if run.log or prot.pdb doesn't exist
            if not os.path.isfile(run_log_file) or not os.path.isfile(prot_pdb_file):
                continue

            # Initialize variables to store the run times
            step1_time = step2_time = step3_time = step4_time = total_time = None
            rot_stat = None

            # Initialize variables for PDB details
            residues = set()
            non_protein_molecules = set()
            water_molecules = set()

            # Read the run.log file
            with open(run_log_file, 'r') as file:
                for line in file:
                    if "Total time for step1" in line:
                        step1_time = round(int(line.split()[-2]) / 60)  # Convert seconds to minutes and round
                    elif "Total time for step2" in line:
                        step2_time = round(int(line.split()[-2]) / 60)
                    elif "Total time for step3" in line:
                        step3_time = round(int(line.split()[-2]) / 60)
                    elif "Total time on MC" in line:
                        step4_time = round(int(line.split()[-2]) / 60)

            # Read the rot_stat file if it exists
            if os.path.isfile(rot_stat_file):
                with open(rot_stat_file, 'r') as file:
                    for line in file:
                        if line.strip().startswith("Total"):
                            parts = line.split()
                            if parts:  # Check if line contains data
                                rot_stat = int(parts[-1])  # Convert the last value to integer
                            break  # Stop reading after finding the Totals line

            # Analyze the prot.pdb file
            with open(prot_pdb_file, 'r') as file:
                for line in file:
                    if line.startswith("ATOM"):  # Protein residues
                        chain_id = line[21].strip()  # Chain ID
                        residue_seq = line[22:26].strip()  # Residue sequence number
                        residue_name = line[17:20].strip()  # Residue name
                        residues.add((chain_id, residue_seq, residue_name))
                    elif line.startswith("HETATM"):
                        chain_id = line[21].strip()  # Chain ID
                        residue_seq = line[22:26].strip()  # Residue sequence number
                        residue_name = line[17:20].strip()  # Residue name
                        if residue_name == "HOH":  # Count water molecules
                            water_molecules.add((chain_id, residue_seq, residue_name))
                        else:  # Count non-protein molecules
                            non_protein_molecules.add((chain_id, residue_seq, residue_name))

            # Count unique residues, non-protein molecules, and water molecules
            num_residues = len(residues)
            num_non_protein_molecules = len(non_protein_molecules)
            num_waters = len(water_molecules)

            # If any times are found, add to data list
            if step1_time is not None and step2_time is not None and step3_time is not None and step4_time is not None:
                total_time = step1_time + step2_time + step3_time + step4_time
                data.append([
                    pdb_folder,
                    num_residues,
                    num_non_protein_molecules,
                    num_waters,
                    rot_stat,
                    step1_time,
                    step2_time,
                    step3_time,
                    step4_time,
                    total_time
                ])

    # Convert the data list to a dataframe
    df = pd.DataFrame(
        data,
        columns=[
            'PDB',
            'Residues',
            'Cofactors',
            'Waters',
            'Conformers',
            'Step1 (premcce)',
            'Step2 (rotamer making)',
            'Step3 (PBE + vdw)',
            'Step4 (Monte Carlo)',
            'Total Time (minutes)'
        ]
    )

    # Sort the dataframe by 'Conformers' column in ascending order
    df = df.sort_values(by='Conformers', ascending=True)

    # Save dataframe to CSV
    df.to_csv('run_times.csv', index=False)
    print(f"Run times saved to 'run_times.csv'.")

    # Save dataframe to TXT in the form of a pandas DataFrame
    with open('run_times.txt', 'w') as f:
        f.write(df.to_string(index=False))
    print(f"Run times saved to 'run_times.txt'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract run times and rot_stat from run.log files.")
    parser.add_argument("runs_folder", nargs="?", default="runs", help="Directory containing pdb folders with run.log and prot.pdb files")
    args = parser.parse_args()

    extract_run_times(args.runs_folder)

