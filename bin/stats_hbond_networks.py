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
        return None, None, None

    # Build an undirected graph from the path
    G = nx.Graph()
    for i in range(len(nodes) - 1):
        G.add_edge(nodes[i], nodes[i + 1])

    # Find the largest connected component
    largest_component = max(nx.connected_components(G), key=len)

    # Compute SHA256 hash from the original network line
    network_hash = hashlib.sha256(network_line.encode()).hexdigest()

    return network_hash, network_line, len(largest_component)

def format_network_with_percentages(network_edges, edge_percentages):
    """Formats network structure with percentages properly aligned below each connection."""
    network_line = " -> ".join(network_edges)  # Network path in a single line
    percentage_line = "PW %: "  # Start with "PW%: "
    total_prob = 1.0
    position_tracker = len(percentage_line)

    for i in range(len(network_edges) - 1):
        edge = f"{network_edges[i]} -> {network_edges[i + 1]}"
        edge_prob = edge_percentages.get(edge, 0) / 100  # Convert to decimal
        total_prob *= edge_prob if edge_prob > 0 else 1  # Avoid zero probability

        # Find the position of '->' in the network line
        arrow_position = network_line.find("->", position_tracker)
        space_padding = arrow_position - len(percentage_line)  # Ensure it's under '->'

        # Append spaces and percentage (now with a % sign at the end)
        percentage_line += " " * space_padding + f"{edge_percentages.get(edge, 0):.2f}%"
        position_tracker = arrow_position + 2  # Move past "->" for next search

    # Format total probability as exponential notation and append the % sign
    uncorr_percentage = f"{total_prob * 100:.2e}%"  # Scientific notation for Uncorr %

    return f"Network structure:\n{network_line}\n{percentage_line}\nUncorr %: {uncorr_percentage}\n"

def analyze_networks(input_dir, top_n, node_min):
    """Analyzes network similarities across all microstate PDBs in the directory."""
    if not os.path.exists(input_dir):
        print(f"Error: Directory '{input_dir}' not found.")
        return

    output_file_name = f"MCCE-ms_hbond_network_stats_nodemin{node_min}.txt"
    output_file_path = os.path.join(input_dir, output_file_name)

    network_counts = defaultdict(int)
    network_examples = {}
    network_files = defaultdict(list)
    edge_counts = defaultdict(int)
    total_microstate_pdbs = 0
    found_valid_networks = False

    for filename in os.listdir(input_dir):
        if "MCCE-ms_hbond_network_stats" in filename or not filename.endswith(".txt"):
            continue

        total_microstate_pdbs += 1
        filepath = os.path.join(input_dir, filename)
        try:
            with open(filepath, "r") as file:
                network_lines = [line.strip() for line in file if line.strip()]

                for network_line in network_lines:
                    if " -> " not in network_line:
                        continue

                    network_hash, network_str, node_count = hash_network(network_line)

                    if network_hash and node_count >= node_min:
                        network_counts[network_hash] += 1
                        network_files[network_hash].append(filename)
                        if network_hash not in network_examples:
                            network_examples[network_hash] = network_str
                        found_valid_networks = True

                        # Count pairwise connections
                        nodes = network_str.split(" -> ")
                        for i in range(len(nodes) - 1):
                            edge = f"{nodes[i]} -> {nodes[i + 1]}"
                            edge_counts[edge] += 1

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

        for i, (network_hash, count) in enumerate(sorted_networks[:top_n]):
            ratio = count / total_microstate_pdbs
            percentage = ratio * 100

            if ratio == 1:
                E = 0
                k = A
            else:
                E = -1.364 * math.log10(ratio)
                k = A * 10**(-E/1.364)

            # Compute edge percentages
            network_edges = network_examples[network_hash].split(" -> ")
            edge_percentages = {edge: (edge_counts[edge] / total_microstate_pdbs) * 100 for edge in edge_counts}

            formatted_network = format_network_with_percentages(network_edges, edge_percentages)

            example_files = ", ".join(network_files[network_hash])
            result = (
                f"Network {i+1}: {count} microstate PDBs ({percentage:.2f}%) share this network.\n"
                f"Network Rate & Energy: k = {k:.2e} /sec, E = {E:.2e} kcal/mol \n"
                f"{formatted_network}\n"
                f"Found in files: {example_files}\n"
                f"{'-' * 200}\n"
            )
            print(result)
            output_file.write(result)

        if sorted_networks:
            most_common_percentage = (sorted_networks[0][1] / total_microstate_pdbs) * 100
            summary = f"Most common network is shared by {most_common_percentage:.2f}% of microstate PDBs."
            print(summary)
            output_file.write(summary)

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
        default=5,
        help="Minimum number of nodes for networks to be considered (default: 5)"
    )

    args = parser.parse_args()
    analyze_networks(args.input_dir, args.topnets, args.node_min)

