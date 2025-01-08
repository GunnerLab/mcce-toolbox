#!/usr/bin/env python
"""
Created on Tues Dec 03 16:01:15 2024

@author: Gehan Ranepura & Koreena Sookhai
"""

import ast
import argparse
import networkx as nx
import os
import re

# Function to process the input data
def process_conserved_graph(file_path):
    with open(file_path, 'r') as file:
        input_text = file.read()

    # Define patterns to handle both "conserved" and non-conserved versions. Find start positions for nodes and edges
    node_pattern = re.compile(r"List of (conserved )?nodes:")
    edge_pattern = re.compile(r"List of (conserved )?edges:")
    nodes_match = node_pattern.search(input_text)
    edges_match = edge_pattern.search(input_text)

    if nodes_match:
        nodes_start = nodes_match.end()
    else:
        raise ValueError("Node list not found in the input text.")

    if edges_match:
        edges_start = edges_match.start()
    else:
        raise ValueError("Edge list not found in the input text.")

    # Extract nodes and convert to list
    nodes_text = input_text[nodes_start:edges_start].strip()
    try:
        nodes = ast.literal_eval(nodes_text.replace("'", '"').replace(" ", ","))
    except Exception as e:
        raise ValueError(f"Error parsing nodes: {e}")
    print(f"Nodes:\n{nodes}\n")

    # Extract edges and convert to list of tuples
    edges_text = input_text[edges_start + len(edges_match.group(0)):].strip()
    try:
        edges_list = ast.literal_eval(edges_text.replace("'", '"').replace(" ", ","))
        edges = [tuple(edge) for edge in edges_list]
    except Exception as e:
        raise ValueError(f"Error parsing edges: {e}")
    print(f"Edges:\n{edges}\n")

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

# Function to save edges to a file
def save_edges_to_file(edges, output_file, input_filename):
    with open(output_file, 'w') as file:
        # Add title at the top of the file
        #file.write(f"Edges for the argument file: {input_filename}\n\n")
        #file.write("Residue1\tResidue2\n")  # Header for the columns
        
        # Write each edge as two columns
        for edge in edges:
            file.write(f"{edge[0]}\t{edge[1]}\n")

# Main function
def main():
    parser = argparse.ArgumentParser(description="Process conserved H-bond graph data.")
    parser.add_argument('filename', help="Input textfile containing the cgraphs H-bond graph data.")
    args = parser.parse_args()

    # Process the input file to extract nodes and edges from the graph
    base_filename = os.path.basename(args.filename).split('.')[0]  # Get the base name of the input file without extension
    nodes, edges = process_conserved_graph(args.filename)
    output_edges_file = f"Edges_{base_filename}.txt"               # Output filename for edges
    save_edges_to_file(edges, output_edges_file, args.filename) 
    print(f"Edges data has been saved to '{output_edges_file}'.")

    # Find all possible paths between residues in the graph
    all_paths = find_all_paths(nodes, edges)
    output_paths_file = f"FullResPaths_{base_filename}.txt"        # Output filename for paths
    save_paths_to_file(all_paths, output_paths_file, args.filename)
    print(f"Paths data has been saved to '{output_paths_file}'.")

    # Optionally, output the paths to the console as well
    print(f"\n===============================================================================================================================")
    print(f"\nFull Paths Between Residues (At Least Two Residues):")
    for start_residue, paths in all_paths.items():
        print(f"\nPaths starting from {start_residue}:")
        for path in paths:
            print(" -> ".join(path))

if __name__ == "__main__":
    main()

