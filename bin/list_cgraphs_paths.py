#!/usr/bin/env python
"""
Created on Tues Dec 03 16:01:15 2024

@author: Gehan Ranepura & Koreena Sookhai
"""

import ast
import argparse
import networkx as nx
import os

# Function to process the input data
def process_conserved_graph(file_path):
    with open(file_path, 'r') as file:
        input_text = file.read()
    
    # Extract nodes and edges from the file content
    nodes_start = input_text.find("List of conserved nodes:") + len("List of conserved nodes:")
    edges_start = input_text.find("List of conserved edges:")
    
    # Extract nodes and convert to list
    nodes_text = input_text[nodes_start:edges_start].strip()
    nodes = ast.literal_eval(nodes_text.replace("'", '"').replace(" ", ","))
    
    # Extract edges and convert to list of tuples
    edges_text = input_text[edges_start + len("List of conserved edges:"):].strip()
    edges_list = ast.literal_eval(edges_text.replace("'", '"').replace(" ", ","))
    edges = [tuple(edge) for edge in edges_list]
    
    return nodes, edges

# Function to find all paths using networkx
def find_all_paths(nodes, edges):
    # Create a graph from edges
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    # Find and output all paths with at least two residues
    all_paths = {}
    for node in nodes:
        paths_for_node = []
        for target in nodes:
            if node != target:  # Exclude paths from a node to itself
                paths = list(nx.all_simple_paths(G, source=node, target=target))
                for path in paths:
                    if len(path) > 1:  # Only include paths with at least two residues
                        paths_for_node.append(path)
        if paths_for_node:  # Only store nodes that have paths to other nodes
            all_paths[node] = paths_for_node
    
    return all_paths

# Function to save the paths to a file
def save_paths_to_file(all_paths, output_file, input_filename):
    with open(output_file, 'w') as file:
        # Add title at the top of the file
        file.write(f"Full paths for the argument file: {input_filename}\n\n")
        
        # Write paths for each starting residue
        for start_residue, paths in all_paths.items():
            file.write(f"Paths starting from {start_residue}:\n")
            for path in paths:
                file.write(" -> ".join(path) + "\n")
            file.write("\n")

# Main function
def main():
    parser = argparse.ArgumentParser(description="Process conserved H-bond graph data.")
    parser.add_argument('filename', help="Input file containing the conserved H-bond graph data.")
    args = parser.parse_args()

    # Process the input file
    nodes, edges = process_conserved_graph(args.filename)
    all_paths = find_all_paths(nodes, edges)
    
    # Create dynamic output file name and save the paths
    base_filename = os.path.basename(args.filename).split('.')[0]  # Get the base name of the input file without extension
    output_file = f"FullResPaths_{base_filename}.txt"  # Create the output filename
    save_paths_to_file(all_paths, output_file, args.filename)

    # Optionally, output the paths to the console as well
    print(f"Full Paths Between Residues (At Least Two Residues):")
    for start_residue, paths in all_paths.items():
        print(f"\nPaths starting from {start_residue}:")
        for path in paths:
            print(" -> ".join(path))

if __name__ == "__main__":
    main()

