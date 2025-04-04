#!/usr/bin/env python
"""
Created on Mar 16 06:00:00 2025

@author: Gehan Ranepura
"""

import os
import glob
import shutil
import re
import networkx as nx
import argparse

# Argument parser for input directory and residue pair file
parser = argparse.ArgumentParser(description="Process hydrogen bond graph networks from run_detect_hbond.py output.")
parser.add_argument("-i",         type=str, default="pdb_output_mc_hbonds", help="Input directory containing hbond text files (default: pdb_output_mc_hbonds)")
parser.add_argument("-resi_list", type=str,                                 help="File containing entry and exit residues of interest (format: XXXCYYYY).")
parser.add_argument("-node_min",  type=int, default=2,                      help="Minimum number of nodes in the graph to process (default: 2)")
args = parser.parse_args()

# Get the input directory and check if it exists
input_dir = args.i
if not os.path.exists(input_dir):
    print(f"Error: Input directory '{input_dir}' does not exist. Please check the path.")
    exit(1)

# Check if output directory exists and prompt the user
OUTPUT_DIR = f"{input_dir}_networks"
if os.path.exists(OUTPUT_DIR):
    user_input = input(f"The output directory '{OUTPUT_DIR}' already exists. Do you want to delete and remake it? (y/n): ").strip().lower()
    if user_input == "y":
        shutil.rmtree(OUTPUT_DIR)  # Delete directory
        os.makedirs(OUTPUT_DIR)    # Remake directory
    else:
        print(f"Please rename or remove '{OUTPUT_DIR}' before running the script again.")
        exit()
else:
    os.makedirs(OUTPUT_DIR)  # Create directory if it doesn't exist

# Function to parse residue information (ENTRY/EXIT residues)
def parse_entryexit_info(residue_info):
    residue_name = residue_info[:3]      # Charcters  1-3 are residue name
    chain_id = residue_info[3]           # Character  4 is the chain ID
    residue_number = residue_info[4:9]   # Characters 5-8 are the residue number
    return residue_name, chain_id, residue_number

# Function to parse residue information (DONOR/ACCEPTOR residues)
def parse_donoracceptor_info(residue_info):
    residue_name = residue_info[:3]      # Characters 1-3 are residue name
    chain_id = residue_info[5]           # Character  6 is the chain ID
    residue_number = residue_info[6:10]  # Characters 7-10 are the residue number
    return residue_name, chain_id, residue_number

# Function to parse ENTRY/EXIT residues in the resi_list
def read_residue_pairs(file_path):
    entry_residues_i = set()
    exit_residues_i = set()
    with open(file_path, 'r') as file:
        lines = file.readlines()

        # Skip the header row (first line)
        for line in lines[1:]:     # Start from the second line
            parts = line.strip().split()
            if len(parts) == 1:    # Entry residue only
                entry_residues_i.add(parts[0])
            elif len(parts) == 2:  # Entry and Exit residue
                entry_residues_i.add(parts[0])
                exit_residues_i.add(parts[1])

    # Parse entry and exit residues
    entry_residues = {parse_entryexit_info(residue) for residue in entry_residues_i}
    exit_residues = {parse_entryexit_info(residue) for residue in exit_residues_i}

    return entry_residues, exit_residues

# Function to process an individual hydrogen bond network file
# Function to process an individual hydrogen bond network file
def process_hbond_graph(file_path, entry_residues, exit_residues):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    nodes = set()
    edges = []

    ENTRY = set()
    EXIT = set()

    # Process each line and build graph edges based on donor/acceptor pairs
    for line in lines:
        parts = line.split()
        if len(parts) < args.node_min:
            continue  # Skip malformed lines

        donor, acceptor = parts[0], parts[1]
        nodes.update([donor, acceptor])  # Create nodes
        edges.append((donor, acceptor))  # Create edges

        donor_residue_info = donor.split('_')[0]  # Extract the part before the conformer info
        acceptor_residue_info = acceptor.split('_')[0]

        donor_residue_name, donor_chain_id, donor_residue_number = parse_donoracceptor_info(donor_residue_info)
        acceptor_residue_name, acceptor_chain_id, acceptor_residue_number = parse_donoracceptor_info(acceptor_residue_info)

        # Match donor and acceptor residues against entry/exit residues
        if (donor_residue_name, donor_chain_id, donor_residue_number) in entry_residues:
            ENTRY.add(donor)     # Add donor to ENTRY list if it's in entry_residues
        if (acceptor_residue_name, acceptor_chain_id, acceptor_residue_number) in entry_residues:
            ENTRY.add(acceptor)  # Add acceptor to ENTRY list if it's in entry_residues

        if (donor_residue_name, donor_chain_id, donor_residue_number) in exit_residues:
            EXIT.add(donor)     # Add donor to EXIT list if it's in exit_residues
        if (acceptor_residue_name, acceptor_chain_id, acceptor_residue_number) in exit_residues:
            EXIT.add(acceptor)  # Add acceptor to EXIT list if it's in exit_residues

    if len(nodes) < args.node_min:
        return None, None, None, None  # Skip this network if it has fewer than 2 nodes

    return list(nodes), list(edges), list(ENTRY), list(EXIT)


# Function for alphauerical path sorting
def natural_sort_key(s):
    """Extracts numeric and non-numeric parts to sort alphanumerically using regex."""
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]

# Function to build graph and save paths between entry and exit residues
def build_graph(ENTRY, EXIT, nodes, edges, output_file):
    # Create the graph object
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    all_paths = set()  # Dictionary to store paths for each node

    # Check if the graph has at least 2 nodes
    if len(G.nodes) < args.node_min:
        print(f"Skipping network with fewer than node_min.")
        return  # Exit the function if the graph has fewer than 2 nodes

    all_paths = set()  # Dictionary to store paths for each node

    # Find all simple paths between specified entry and exit residues
    for entry_resi in ENTRY:
        for exit_resi in EXIT:
            print(f"Entry/Exit: {entry_resi} <--> {exit_resi}")
            # Find all simple paths starting from entry_residue and ending at exit_residue
            try:
                paths = list(nx.all_simple_paths(G, source=entry_resi, target=exit_resi))
                for path in paths:
                    if len(path) > 1:  # Only add paths that have at least two nodes
                       all_paths.add(tuple(path))  # Convert list to tuple for uniqueness
            except nx.NetworkXNoPath:
               continue  # Skip if no path exists

    # Find all simple paths between entry resi GLU-1P0041_005 and exit resi GLU-1P0050_005
    #GLU-1P0041_005 -> ARG+1P0048_006 -> GLU-1P0050_005
    #try:
    #    paths = list(nx.all_simple_paths(G, source="GLU-1P0041_005", target="GLU-1P0050_005"))
    #    for path in paths:
    #         all_paths.add(tuple(path))  # Convert list to tuple for uniqueness
    #     except nx.NetworkXNoPath:
    #            continue  # Skip if no path exists

    # Find all simple paths leading to specified exit residues
    #for node in nodes:
    #    for exit_resi in EXIT:
    #        print(f"Exit: {exit_resi}")
    #        try:
    #            paths = list(nx.all_simple_paths(G, source=node, target=exit_resi))
    #            for path in paths:
    #                all_paths.add(tuple(path))  # Convert list to tuple for uniqueness
    #        except nx.NetworkXNoPath:
    #               continue  # Skip if no path exists

    # Save unique and alphanumercally sorted paths to output file
    sorted_paths = sorted(all_paths, key=lambda path: [natural_sort_key(node) for node in path])
    try:
        with open(output_file, 'w') as file:
            for path in sorted_paths:
                file.write(" -> ".join(path) + "\n")
        print(f"Unique sorted paths saved to {output_file}")
    except Exception as e:
        print(f"Error writing to file {output_file}: {e}")


# Main function to process all files and build graphs
def process_directory(directory, entry_residues, exit_residues):
    files = glob.glob(os.path.join(directory, "*.txt"))  # Only get .txt files
    files = [f for f in files if not f.endswith("blocking.txt")]

    if not files:
        print(f"No valid hydrogen bond files found in '{directory}'.")
        return

    for file_path in files:
        base_filename = os.path.basename(file_path).split('.')[0]
        print(f"Processing file {base_filename} ...")

        # Extract nodes and edges from the hydrogen bond graph file
        nodes, edges, ENTRY, EXIT = process_hbond_graph(file_path, entry_residues, exit_residues)
        print(f"NODES: {nodes}")
        print()
        print(f"EDGES: {edges}")
        print()
   
        print(f"ENTRY: {ENTRY}")
        print(f"EXIT:  {EXIT}")

        # Build and save the graph for this file
        output_file = os.path.join(OUTPUT_DIR, f"Paths_{base_filename}_Resi-EntryToExit.txt")
        build_graph(ENTRY, EXIT, nodes, edges, output_file)
        print("---------------------------------------------------------------------------------------------------------")

# Ensure residue pairs file is provided and exists
if not args.resi_list or not os.path.exists(args.resi_list):
    print("Error: Residue pair list file is missing or does not exist.")
    exit(1)

# Read the residue pairs from the file
entry_residues, exit_residues = read_residue_pairs(args.resi_list)

# Run the script
process_directory(input_dir, entry_residues, exit_residues)

