#!/usr/bin/env python
"""
Created on Feb 11 01:25:00 2025

@author: Gehan Ranepura
"""

import numpy as np
import math
import random
import sys
import argparse

# This python script creates N-number of water conformers as the output "HOHs_confs.pdb" in the format of MCCE step2_out.pdb
# An argument of a PDB file containing an oxygen atom is used to add two hydrogen atoms to create a water molecule
# The two hydrogen atoms make a specified bond length and bond angle with the oxygen atom

# Define command-line arguments
parser = argparse.ArgumentParser(description="Generate water conformers for an oxygen atom of PDB file.")
parser.add_argument("input_pdb", type=str, help="The input step1_out.pdb file containing the oxygen atom.")
parser.add_argument("-N", type=int, default=25, help="Number of water conformers to generate (default: 25).")
args = parser.parse_args()

# Read the input PDB file to extract the oxygen atom coordinates and the ChainRes_ID
with open(args.input_pdb) as file:
    for line in file:
        if line.startswith("HETATM") and line[17:20].strip() == "HOH" and line[12:16].strip() == "O":
            newline = line.strip().split()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            oxygen_coordinates = [x, y, z]

            Resi    = str(line[17:20])
            ChainID = str(line[21])
            ResID   = str(line[22:26])

            #f = open(args.input_pdb)
            #t = f.readline().split()
            #b1 = str(t[4])
            #b2 = float(t[5])
            #ChainRes_ID = b1 + "%04d" % (b2,) + "_"

            ChainRes_ID = ChainID + ResID + "_"
            step1_water = Resi + " " + ChainID + ResID
            print("\n") 
            print(step1_water, "---> Oxygen Coor  = ", [x, y, z])
            #print("O_coord  = ", [x, y, z])


# Calculate the hydrogen positions for N conformers of the water molecule
# bond_length =   0.96 # The typical bond length for an oxygen-hydrogen bond in a water molecule
# bond_angle  = 104.5  # The typical bond angle for a water molecule

with open("HOH_confs.pdb", 'w') as file:
    for i in range(args.N):
        def add_Hs(O_coord, bond_length=0.96, bond_angle=104.5):
            """Add two hydrogen atoms to a PDB file with an oxygen atom."""

            # Add first hydrogen atom randomly within a spherical radius of 1 angstrom
            rand_dir = np.random.rand(3)
            rand_dir /= np.linalg.norm(rand_dir)
            z = [np.random.choice([-1, 1]), np.random.choice([-1, 1]), np.random.choice([-1, 1])]
            H1_coord = O_coord + bond_length * rand_dir * z

            # Calculate unit vector along O-H1 bond
            bond_vector = H1_coord - O_coord
            bond_unit_vector = bond_vector / np.linalg.norm(bond_vector)

            # Calculate an orthogonal vector to O-H1 bond
            orthogonal_vector = np.cross(bond_vector, [0, 0, 1])
            if np.allclose(orthogonal_vector, 0):
                orthogonal_vector = np.cross(bond_vector, [0, 1, 0])
            orthogonal_unit_vector = orthogonal_vector / np.linalg.norm(orthogonal_vector)

            # Calculate second hydrogen atom position
            bond_angle = bond_angle * math.pi / 180
            a = np.cos(bond_angle / 2)
            b, c, d = orthogonal_unit_vector * np.sin(bond_angle / 2)
            q = np.array([a, b, c, d])
            R = np.array([[a ** 2 + b ** 2 - c ** 2 - d ** 2, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                          [2 * (b * c + a * d), a ** 2 + c ** 2 - b ** 2 - d ** 2, 2 * (c * d - a * b)],
                          [2 * (b * d - a * c), 2 * (c * d + a * b), a ** 2 + d ** 2 - b ** 2 - c ** 2]])
            H2_coord = O_coord + np.matmul(R, bond_unit_vector) * bond_length

            # Define a rotation matrix around the O atom
            theta = np.random.uniform(0, 2 * math.pi)
            phi = np.random.uniform(0, math.pi)

            R_x = np.array([[1, 0, 0],
                            [0, math.cos(theta), -math.sin(theta)],
                            [0, math.sin(theta), math.cos(theta)]])

            R_y = np.array([[math.cos(phi), 0, math.sin(phi)],
                            [0, 1, 0],
                            [-math.sin(phi), 0, math.cos(phi)]])

            Rot = np.dot(R_y, R_x)

            # Randomly rotate the two hydrogen positions
            p1 = np.array(O_coord)
            p2 = np.array(H1_coord)
            p3 = np.array(H2_coord)

            H1_coord = np.dot(Rot, p2 - p1) + p1
            H2_coord = np.dot(Rot, p3 - p1) + p1

            # Print the Bond Angle between HOH, the Bond Lengths the between O and H atoms, and the coordinates of H1 and H2
            v1 = np.array(H1_coord) - np.array(O_coord)
            v2 = np.array(H2_coord) - np.array(O_coord)
            HOHangle = np.degrees(np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))
            h1_L = np.linalg.norm(np.array(O_coord) - np.array(H1_coord))
            h2_L = np.linalg.norm(np.array(O_coord) - np.array(H2_coord))
            H1 = ['{:07.3f}'.format(H1_coord[0]), '{:07.3f}'.format(H1_coord[1]), '{:07.3f}'.format(H1_coord[2])]
            H2 = ['{:07.3f}'.format(H2_coord[0]), '{:07.3f}'.format(H2_coord[1]), '{:07.3f}'.format(H2_coord[2])]

            print("HOH_{} --> Angle = {}, O-H1 = {}, O-H2 = {}, H1 = {}, H2 = {}".format('{:03d}'.format(i+1), '{:.1f}'.format(HOHangle),
                                                                                         '{:.2f}'.format(h1_L), '{:.2f}'.format(h2_L),
                                                                                          H1, H2))

            return H1_coord, H2_coord


        O_coord = [oxygen_coordinates[0], oxygen_coordinates[1], oxygen_coordinates[2]]
        H1_coord, H2_coord = add_Hs(O_coord)
        h1_x, h1_y, h1_z = H1_coord[0], H1_coord[1], H1_coord[2]
        h2_x, h2_y, h2_z = H2_coord[0], H2_coord[1], H2_coord[2]


        # Write Final Coordinates
        O  = O_coord
        H1 = [h1_x, h1_y, h1_z]
        H2 = [h2_x, h2_y, h2_z]

        # Write pdb
        #AtomO  = '{:07.3f}'.format(O[0])+""+'{:07.3f}'.format(O[1])+""+'{:07.3f}'.format(O[2])
        #AtomH1 = '{:07.3f}'.format(H1[0])+""+'{:07.3f}'.format(H1[1])+""+'{:07.3f}'.format(H1[2])
        #AtomH2 = '{:07.3f}'.format(H2[0])+""+'{:07.3f}'.format(H2[1])+""+'{:07.3f}'.format(H2[2])
        AtomO  = "{:8.3f}{:8.3f}{:8.3f}".format(O[0], O[1], O[2])
        AtomH1 = "{:8.3f}{:8.3f}{:8.3f}".format(H1[0], H1[1], H1[2])
        AtomH2 = "{:8.3f}{:8.3f}{:8.3f}".format(H2[0], H2[1], H2[2])

        #file.write("HETATM    1  O   HOH "+ChainRes_ID+""+'{:03}'.format(i+1)+" "+AtomO+ "  1.520      -0.800      01O000M000"+"\n")
        #file.write("HETATM    2  H1  HOH "+ChainRes_ID+""+'{:03}'.format(i+1)+" "+AtomH1+"  1.100       0.400      01O000M000"+"\n")
        #file.write("HETATM    3  H2  HOH "+ChainRes_ID+""+'{:03}'.format(i+1)+" "+AtomH2+"  1.100       0.400      01O000M000"+"\n")
        file.write("HETATM    1  O   HOH {}{:03}{}".format(ChainRes_ID, i+1, AtomO)  + "   1.520      -0.800      01O000M000\n")
        file.write("HETATM    2  H1  HOH {}{:03}{}".format(ChainRes_ID, i+1, AtomH1) + "   1.100       0.400      01O000M000\n")
        file.write("HETATM    3  H2  HOH {}{:03}{}".format(ChainRes_ID, i+1, AtomH2) + "   1.100       0.400      01O000M000\n")

 

