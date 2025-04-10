#!/usr/bin/env python
"""
Created on Apr 09 05:20:00 2025

@author: Gehan Ranepura
"""

import os
import argparse
import math
from collections import Counter, defaultdict

def parse_network_file(filepath, node_min):
    """Extract all valid networks (as node lists) from a given file."""
    networks = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or "->" not in line:
                continue
            nodes = [n.strip() for n in line.split("->")]
            if len(nodes) >= node_min:
                networks.append(nodes)
    return networks

def print_and_write(output_file, text):
    """Print to console and write to file."""
    print(text)
    output_file.write(text + '\n')

def main(input_dir, topnets, node_min, A, output_file):
    all_networks = []
    all_pairwise_presence = defaultdict(set)
    network_file_map = defaultdict(set)

    filenames = [f for f in os.listdir(input_dir)
                 if f.startswith("Paths_ms_pdb_") and f.endswith("_Resi-EntryToExit.txt")]
    total_files = len(filenames)

    with open(output_file, 'w') as out_file:
        print_and_write(out_file, f"\nTotal microstate PDBs analyzed: {total_files}\n")
        print_and_write(out_file, f"Top {topnets} networks (with at least {node_min} nodes):\n")

        file_networks_dict = {}

        for filename in filenames:
            full_path = os.path.join(input_dir, filename)
            file_networks = parse_network_file(full_path, node_min)
            file_networks_strs = []

            for path_nodes in file_networks:
                network_str = " -> ".join(path_nodes)
                file_networks_strs.append(network_str)
                all_networks.append(network_str)
                network_file_map[network_str].add(filename)

                for i in range(len(path_nodes) - 1):
                    pair = (path_nodes[i], path_nodes[i + 1])
                    all_pairwise_presence[pair].add(filename)

            file_networks_dict[filename] = file_networks_strs

        network_counts = Counter(all_networks)
        top_networks = network_counts.most_common(topnets)

        for idx, (network, count) in enumerate(top_networks, start=1):
            nodes = network.split(" -> ")
            arrows = ["->"] * (len(nodes) - 1)

            # Build network string and collect arrow positions
            parts = [nodes[0]]
            arrow_starts = []
            for i in range(len(arrows)):
                parts.append(f" {arrows[i]} {nodes[i + 1]}")
                arrow_start_pos = sum(len(p) for p in parts[:-1]) + 1
                arrow_starts.append(arrow_start_pos)
            network_str = ''.join(parts)

            network_percentage = (count / total_files) * 100
            ratio = count / total_files
            if ratio == 1:
                E = 0
                k = A
            else:
                E = -1.364 * math.log10(ratio)
                k = A * 10 ** (-E / 1.364)

            print_and_write(out_file, f"{idx}. Count: {count} microstate PDBs ({network_percentage:6.2f}%) have this fully connected network.")
            print_and_write(out_file, f"Network Rate & Energy: k = {k:.2e} /sec, E = {E:.2e} kcal/mol")
            print_and_write(out_file, f"Network: {network_str}")

            # Build PW % line
            pw_percents = []
            for i in range(len(nodes) - 1):
                pair = (nodes[i], nodes[i + 1])
                percent = len(all_pairwise_presence[pair]) / total_files * 100
                pw_percents.append(f"{percent:6.2f}%")

            pw_line = [' '] * len(network_str)
            for start, percent in zip(arrow_starts, pw_percents):
                for j, ch in enumerate(percent):
                    if 0 <= start + j < len(pw_line):
                        pw_line[start + j - 1] = ch
            print_and_write(out_file, "PW %:    " + ''.join(pw_line))

            # Build Subnet % line
            subnet_percents = []
            for i in range(1, len(nodes)):
                subnet_nodes = nodes[:i+1]
                subnet_str = " -> ".join(subnet_nodes)
                count_subnet = sum(
                    any(subnet_str in net for net in nets)
                    for nets in file_networks_dict.values()
                )
                percent = count_subnet / total_files * 100
                subnet_percents.append(f"{percent:6.2f}%")

            subnet_line = [' '] * len(network_str)
            for start, percent in zip(arrow_starts, subnet_percents):
                for j, ch in enumerate(percent):
                    if 0 <= start + j < len(subnet_line):
                        subnet_line[start + j - 1] = ch
            print_and_write(out_file, "Subnet %:" + ''.join(subnet_line))

            # File list
            print_and_write(out_file, f"\nFound in files: {', '.join(sorted(network_file_map[network]))}")
            print_and_write(out_file, f"{'-' * 200}\n")

    print(f"Results saved to: {output_file}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze network similarities in microstate PDB-associated files.")
    parser.add_argument("-input_dir", type=str,   default="ms_pdb_output_hbonds_networks", help="Directory containing network files (default: %(default)s)")
    parser.add_argument("-topnets",   type=int,   default=5,                               help="Number of top networks to display (default: %(default)s)")
    parser.add_argument("-node_min",  type=int,   default=5,                               help="Minimum number of nodes for networks to be considered (default: %(default)s)")
    parser.add_argument("-A",         type=float, default=10**13,                          help="Pre-exponential factor (default: %(default)s)")
    args = parser.parse_args()

    output_file_name = f"MCCE-ms_hbond_network_stats_nodemin{args.node_min}.txt"
    output_file_path = os.path.join(args.input_dir, output_file_name)

    main(args.input_dir, args.topnets, args.node_min, args.A, output_file_path)

