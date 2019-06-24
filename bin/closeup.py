#!/usr/bin/env python
"""
This program gives a list of atoms that are within 4 angstroms of an OMC
"""

import sys
import numpy as np

# read pdb file
Six_GIX_lines = open(sys.argv[1]).readlines()
pdb_name = Six_GIX_lines[0][62:66] # get the 4 letter code of protein from pdb file
pdb_l = [pdb_name]   # put the 4 letter code name in a list
b = ['.pdb'] #create list with file name
pdb = pdb_l + b # concatenate the list
pdb_out = ''.join(pdb) # join the two lists into a string


 

# make a list with all atoms: ATOM or HETATM
atm_lst = [lines for lines in Six_GIX_lines if lines[:4]== "ATOM" or lines[:6] == "HETATM" ]

#make a list with all OMC atoms
OMC = [lines for lines in atm_lst if lines[17:20] == "OMC"]

#make a list with all OMC atoms
trg_atm = "OMC"
OMC = [atoms for atoms in atm_lst if atoms[13:16] == "OMC"]

#extract coordinates of OMC
coor_OMC = [lines[31:54] for lines in OMC]

# convert string elemens in above list to floats--gives list with floats, not grouped by ATOM
pre_loc_OMC = [float(j) for i in coor_OMC for j in i.split()]

#make list with elements that are grouped by atoms position
pos_OMC = [pre_loc_OMC[i:i+3] for i in range(0,len(pre_loc_OMC),3)]



#repeat the above procedure with all other atoms: make a list of all other atoms'names
#and make a list with their respectivce coordinates

#make a list with all other atoms
other_atms = [atoms for atoms in atm_lst if atoms[13:16] != "OMC" and atoms[0:4] != "TER"]

#extract coordinates of all other atoms
coor_othrs = [lines[31:54] for lines in other_atms]

pre_loc_othrs = [float(j) for i in coor_othrs for j in i.split()] #converting to float, not grouped

pos_othrs = [pre_loc_othrs[j:j+3] for j in range(0,len(pre_loc_othrs),3)] #group other atom coordinates



#use numpy to subtract coordinates
Othrs_coor = np.array(pos_othrs) #matrix with coordinate of other atoms as rows
OMC_coor0 = np.array(pos_OMC[3]) #make a row matrix of coordinates of one OMC ###########------->>>

#identify what is your target chlorophyll
trg_yll = OMC[3][17:26]#######--------------->>>>>

rel = [o_coor - OMC_coor0 for o_coor in Othrs_coor] #subtract coor of CMC from Othrs_coor
rel_sq = np.square(rel) # square the relative distance 

sum_rel_sq = [np.sum(coor) for coor in rel_sq] #sum the coordinates squared
dist_pass_in = [index for index, value in enumerate(sum_rel_sq) if value <= 16] #list with the indices of those atoms that meet criteria
dist_pass = [element for element in sum_rel_sq if element <= 16] #list with distances that pass

#extract atoms names from other_atoms if they share the same index as those elements in dist_pass_in above
atom_list = [atoms[13:16] for i, atoms in enumerate(other_atms) for elements in dist_pass_in if i== elements]
a_l = [at_nme.strip(' ') for at_nme in atom_list] # strip the list with atom names

#make a list with all information of atoms that meet distance criteria of < 16 ang squared
l = [atoms for i, atoms in enumerate(other_atms) for elements in dist_pass_in if i == elements]

l_p = [rows for rows in l if rows[17:26] != trg_yll] # remove atoms that are part of target chlorophyll


#list with distances associated with non-target chlorophyll passing atoms
coor_non_chlro = [lines[31:54] for lines in l_p]# extract coor of non-'yll passing atoms
p_loc_non_chlro = [float(j) for i in coor_non_chlro for j in i.split()] #converting to float, not grouped
pos_non_chlro = [p_loc_non_chlro[j:j+3] for j in range(0,len(p_loc_non_chlro),3)] #group their atom coordinates

non_chlro_coor = np.array(pos_non_chlro) #matrix with coordinate of other not target atoms
rel_non_chlro = [pos_non_chrlo - OMC_coor0 for pos_non_chrlo in non_chlro_coor] #subtract coor of CMC from Othrs_coor
rel_sq_non_chlro = np.square(rel_non_chlro) # square the relative distance
sum_rel_sq_non_trgt = [np.sum(coor) for coor in rel_sq_non_chlro] #sum the coordinates squared


#make a list for the residue names of the atoms
pass_res_name = [residues[17:26] for residues in l_p]

# make a list for the chain ID's
chain_ids = [chains[21] for i, chains in enumerate(other_atms) for elements in dist_pass_in if i == elements]

#make a list with residue sequence number
res_seq = [res_id[22:26] for i, res_id in enumerate(other_atms) for elements in dist_pass_in if i == elements]
#res_seq_id = [resid.strip(' ') for resid in res_seq]# strip the list with residue sequence id's
r_s_id = [int(i) for i in res_seq] # convert string elements 

#list l below should have the length of purged list--is this a finger in my hand?
#make a lists for output puroses of PDB file name, target atoms, target chlorophyll
l_pdb = [pdb_out for i in l_p]
l_trg_at = [trg_atm for i in l_p]
l_trg_yll = [trg_yll for i in l_p]




#for the target chlrophyll include chain id and residue seq; atoms in Range refer to those that are less
#then 16 Angstroms squared from target atom
print('PDB file       Target Atom      Target Chlorophyll    Atoms in Range    Residue     Distance ang. sqd')
print('-----------------------------------------------------------------------------------------------------')
for i, j, k, z, m, n  in zip(l_pdb, l_trg_at, l_trg_yll , atom_list, pass_res_name, sum_rel_sq_non_trgt):
    print('%s %10s %22s %18s %18s %18.4s'% (i, j, k, z, m,n))
print('\n')
 
 
 
