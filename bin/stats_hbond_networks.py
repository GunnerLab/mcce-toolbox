#!/usr/bin/env python
"""
Created on Mar 16 06:00:00 2025

@author: Gehan Ranepura
"""

import os
import hashlib
import argparse
import math
import networkx as nx
from collections import defaultdict

# Define Arrhenius Rate Equation [ k = A·e**(-E / R·T) = 1.364·probability]:
A = 10**13  # sec^(-1)

def hash_network(network_line):
    """Creates a hash for a network path to uniquely identify it while preserving the original input."""
    
    nodes = network_line.split(" -> ")

    # Skip invalid input
    if len(nodes) <= 1:
        return None, None, None  # Not a valid network

    # Build an undirected graph from the path
    G = nx.Graph()
    for i in range(len(nodes) - 1):
        G.add_edge(nodes[i], nodes[i + 1])

    # Find the largest connected component
    largest_component = max(nx.connected_components(G), key=len)
    
    # Ensure the original network_line is preserved
    connected_network_str = network_line  # Keep it exactly as input
    
    # Compute SHA256 hash from the original network line
    network_hash = hashlib.sha256(network_line.encode()).hexdigest()

    return network_hash, connected_network_str, len(largest_component)


def analyze_networks(input_dir, top_n, node_min):
    """Analyzes network similarities across all microstate PDBs in the directory."""
    if not os.path.exists(input_dir):
        print(f"Error: Directory '{input_dir}' not found.")
        return

    output_file_name = f"MCCE-ms_hbond_network_stats_nodemin{node_min}.txt"
    output_file_path = os.path.join(input_dir, output_file_name)

    network_counts = defaultdict(int)
    network_examples = {}  # Store an example network for each hash
    network_files = defaultdict(list)  # Store files where each network is found
    total_microstate_pdbs = 0
    found_valid_networks = False

    for filename in os.listdir(input_dir):
        # Skip files containing "MCCE-ms_hbond_network_stats" in their name
        if "MCCE-ms_hbond_network_stats" in filename:
            continue

        if filename.endswith(".txt"):  # Process only .txt files
            total_microstate_pdbs += 1
            filepath = os.path.join(input_dir, filename)
            try:
                with open(filepath, "r") as file:
                    network_lines = [line.strip() for line in file if line.strip()]

                    # Process each network line separately
                    for network_line in network_lines:
                        # Only process lines that contain " -> " (i.e., networks with at least two nodes)
                        if " -> " not in network_line:
                            continue

                        # Process valid networks with more than one node
                        network_hash, network_str, node_count = hash_network(network_line)

                        if network_hash and node_count >= node_min:  # Only consider valid networks with the required number of nodes
                            network_counts[network_hash] += 1
                            network_files[network_hash].append(filename)  # Store the file where the network is found
                            if network_hash not in network_examples:
                                network_examples[network_hash] = network_str  # Store an example
                            found_valid_networks = True  # We found at least one valid network

            except Exception as e:
                print(f"Error reading file {filename}: {e}")

    if total_microstate_pdbs == 0:
        print("No valid network files found.")
        return

    if not found_valid_networks:
        print(f"Can't find networks with the node_min limit of {node_min}.")
        return

    sorted_networks = sorted(network_counts.items(), key=lambda x: x[1], reverse=True)

    with open(output_file_path, "w") as output_file:
        output_file.write(f"Total microstate PDBs analyzed: {total_microstate_pdbs}\n\n")
        print(f"Total microstate PDBs analyzed: {total_microstate_pdbs}\n")

        for i, (network_hash, count) in enumerate(sorted_networks[:top_n]):  # Show top N networks
            ratio      = (count / total_microstate_pdbs)
            percentage = ratio * 100
 
            # Return the ratio in the form of 10^x (or fraction x 10^x)
            if ratio == 1:
                 E = 0
                 k = A
            else:
                 E = -1.364 * math.log10(ratio)
                 k = A * 10**(-E/1.364)
 
            example_files = ", ".join(network_files[network_hash])  # Get example files
            result = (
                f"Network {i+1}: {count} microstate PDBs ({percentage:.2f}%) share this network.\n"
                f"Network Rate & Energy: k = {k:.2e} /sec, E = {E:.2e} kcal/mol \n"
                f"Network structure: \n{network_examples[network_hash]} \n\n"
                f"Found in files: {example_files}\n"
                f"{'-' * 200}\n"
            )
            print(result)  # Print to console
            output_file.write(result)  # Save to file

        if sorted_networks:
            most_common_percentage = (sorted_networks[0][1] / total_microstate_pdbs) * 100
            summary = f"Most common network is shared by {most_common_percentage:.2f}% of microstate PDBs."
            print(summary)  # Print to console
            output_file.write(summary)  # Save to file

    print(f"Results saved to: {output_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze network similarities in microstate PDB-associated files.")
    parser.add_argument(
        "input_dir",
        nargs="?",
        default="pdb_output_mc_hbonds_networks",
        help="Directory containing network files (default: pdb_output_mc_hbonds_networks)"
    )
    parser.add_argument(
        "-topnets",
        type=int,
        default=5,
        help="Number of top networks to display (default: 5)"
    )
    parser.add_argument(
        "-node_min",
        type=int,
        default=2,
        help="Minimum number of nodes for networks to be considered (default: 2)"
    )

    args = parser.parse_args()
    analyze_networks(args.input_dir, args.topnets, args.node_min)

