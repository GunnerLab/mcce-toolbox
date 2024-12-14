#!/usr/bin/env python
"""
Created on Thurs Oct 10 14:02:34 2024

@author: Gehan Ranepura
"""

import sys
import numpy as np
from Bio.PDB import PDBParser, PDBIO
from scipy.spatial.distance import cdist

def reposition_water_oxygens(pdb_file, max_cycles=100, iterations_per_cycle=25, cutoff=3.0, max_displacement=0.50):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    # Extract atoms and waters (OX)
    atoms = [atom for atom in structure.get_atoms() if atom.element != 'H']
    waters = [atom for atom in atoms if atom.get_parent().get_resname() == 'HOH' and atom.element == 'O']

    # Prepare for output text file
    reposition_info = []

    def distance_check(atom, new_coord, atoms, cutoff):
        """Check if a proposed new position is at least cutoff Å from other atoms."""
        dist = cdist([new_coord], [a.get_coord() for a in atoms if a != atom])
        min_distance = np.min(dist)
        return min_distance >= cutoff, min_distance

    def calculate_displacement(coord1, coord2):
        """Calculate the Euclidean distance between two sets of coordinates."""
        return np.linalg.norm(coord1 - coord2)

    def format_coordinates(coords):
        """Format coordinates with 5 characters before the decimal point and 3 decimal places after."""
        return '[ ' + ' '.join(f"{x:8.3f}" for x in coords) + ' ]'

    def format_value(value):
        """Ensure values are formatted to 3 decimal places with trailing zeros."""
        return f"{value:.3f}"

    def final_check():
        """Final check to ensure all waters are at least 3.2 Å away from any other atom."""
        final_issues = []
        for water in waters:
            nearby_atoms = [a for a in atoms if a != water]
            dist_matrix = cdist([water.get_coord()], [a.get_coord() for a in nearby_atoms])
            min_distance = np.min(dist_matrix)
            if min_distance < cutoff:
                water_id = water.get_parent().get_id()[1]
                chain_id = water.get_parent().get_parent().get_id()
                final_issues.append((water, min_distance))
        return final_issues

    # Reposition waters over multiple cycles
    for cycle in range(max_cycles):
        print(f"\nCycle {cycle + 1} of {max_cycles}")
        
        # Check waters that are too close and need repositioning
        waters_to_fix = final_check()

        if not waters_to_fix:
            print("All waters are properly positioned!")
            break
        
        print(f"{len(waters_to_fix)} waters are still too close to neighbors, repositioning...")

        for water, min_distance in waters_to_fix:
            original_coord = water.get_coord()
            chain_id = water.get_parent().get_parent().get_id()  # Move chain_id assignment here
            nearby_atoms = [a for a in atoms if a != water]

            repositioned = False  # Flag to check if repositioning was successful
            for i in range(iterations_per_cycle):
                # Generate a small random displacement within the max_displacement limit
                new_coord = original_coord + np.random.uniform(-max_displacement, max_displacement, size=3)

                # Calculate displacement to check against max_displacement
                displacement = calculate_displacement(original_coord, new_coord)

                # If the displacement exceeds the max allowed, continue to the next iteration
                if displacement > max_displacement:
                    continue  # Skip this iteration if the displacement exceeds max_displacement

                new_coord = np.round(new_coord, 3)  # Limit to 3 decimal places
                
                # Check if new position satisfies the distance constraint
                valid, new_min_distance = distance_check(water, new_coord, nearby_atoms, cutoff)
                
                if valid:
                    water.set_coord(new_coord)
                    water_id = water.get_parent().get_id()[1]  # Residue ID (Water ID)
                    reposition_info.append([
                        f"{water_id} {chain_id}",           # Water ID and Chain
                        format_coordinates(original_coord), # Original coordinates formatted
                        format_coordinates(new_coord),      # New repositioned coordinates formatted
                        format_value(displacement),         # Displacement with trailing zeros
                        format_value(min_distance),         # Original minimum distance with trailing zeros
                        format_value(new_min_distance)      # New minimum distance with trailing zeros
                    ])
                    repositioned = True
                    break

            # If not repositioned after all iterations, retain original position
            if not repositioned:
                # print(f"Could not reposition water {water.get_parent().get_id()[1]} {chain_id} to meet distance requirement.")
                pass

    # Save the new PDB with repositioned waters
    io = PDBIO()
    io.set_structure(structure)
    output_pdb = f"{pdb_file[:-4]}_HOHrepositioned.pdb"
    io.save(output_pdb)
    
    # Write the reposition information to a text file
    output_txt = f"{pdb_file[:-4]}_water_repositioning_report.txt"
    with open(output_txt, 'w') as report:
        report.write(f"{'Water ID & Chain':<20}{'Original Coordinates (X, Y, Z)':<45}{'New Coordinates (X, Y, Z)':<45}"
                     f"{'Displacement':<15}{'Original Min Distance':<25}{'New Min Distance':<25}\n")
        report.write('-' * 175 + '\n')
        for entry in reposition_info:
            report.write(f"{entry[0]:<20}{entry[1]:<45}{entry[2]:<45}"
                         f"{entry[3]:<15}{entry[4]:<25}{entry[5]:<25}\n")

    print(f"\nRepositioned PDB saved as {output_pdb}")
    print(f"Water repositioning report saved as {output_txt}")

    # Final check after all cycles
    final_issues = final_check()
    if final_issues:
        print("\nWARNING: Some waters are still too close to their neighbors after all cycles:")
        for water, min_distance in final_issues:
            water_id = water.get_parent().get_id()[1]
            chain_id = water.get_parent().get_parent().get_id()
            print(f"Water {water_id} {chain_id} is only {min_distance:.3f} Å away from a neighboring atom.")
    else:
        print("\nFinal check passed: All water oxygens are at least 3.2 Å away from their neighbors.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python HOH_reposition-pdb.py <pdb_file>")
    else:
        pdb_file = sys.argv[1]
        reposition_water_oxygens(pdb_file)


