import sys

try:
    input_file = sys.argv[1] # take the first terminal argument after the Python script is called
except IndexError:
    print("\nPlease include the PDB file after the executable, e.g. 'ProtInfo_J.py 4lzt.pdb'.\n")
    raise

def line_counter(file_data):  # accept file data and output water, ligand, and amino acid count

    # why not just call the file name directly?

    waters = 0
    ligand = 0
    aa = 0

    for this_line in file_data:

        if "HOH" in this_line and ("HETATM" in this_line or "ANISOU" in this_line):
            waters += 1
        if "HETATM" in this_line and "HOH" not in this_line: # is this how ligands should be found for PDB file? MCCE seems to add ligands
            ligand += 1
        #if any(aa in this_line for aa in standard_aa):  # should return similarly to " grep ^ATOM *.pdb | grep ' CA ' | wc -l", as an example
        #    aa += 1
        if " CA " in this_line:
            aa += 1

    return waters, ligand, aa

step1_output = 'step1_out.pdb'
head1_list = 'head1.lst'
run_log = 'run1.log'

with open(input_file, 'r') as file: # for how many waters, ligands, aa's, were in protein file
    data = file.read()

with open(step1_output, 'r') as file:  # for how many waters, ligands, aa's, post-processing
    mcce_data = file.read()

with open(run_log, 'r') as file: # helps identify where changes have occured
    log_data = file.read()

with open(head1_list, 'r') as file: # not sure how this is useful yet
    head_data = file.read()

data = data.splitlines()
mcce_data = mcce_data.splitlines()
log_data = log_data.splitlines()
head_data = head_data.splitlines()

standard_aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", 
               "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

HOH_count, ligand_count, aa_count = line_counter(data)
HOH_mcce_count, ligand_mcce_count, aa_mcce_count = line_counter(mcce_data)

print("\nNew ProtInfo by Jared! It's still in progress!\n")

print("For the input file " + input_file + ", we find the following:\n")

print("               PDB CNT        MCCE CNT        Difference")
print("Amino Acids    " + str(aa_count) + "            " + str(aa_mcce_count) + "               " + str(abs(aa_count - aa_mcce_count)))
print("Ligands        " + str(ligand_count) + "              " + str(ligand_mcce_count) + "                 " + str(abs(ligand_mcce_count - ligand_count)))  # HETATM that are not waters
print("Waters         " + str(HOH_count) + "             " + str(HOH_mcce_count) + "                 " + str(abs(HOH_count - HOH_mcce_count)))  

# is the protein file multi-chain? If so, do for each chain.

print("\nEditable file that controls ligand stripping:") # appears to be 00always_needed.tpl, so give the pathway to it

for line in log_data:

    if "Labeling" in line and "NTR" in line:
        NTR_line = line
    if "Labeling" in line and "CTR" in line:
        CTR_line = line

print("\nThese residues have been modified:")
print("\nTERMINI:\n")

print(NTR_line)
print(CTR_line)


print("\nLIGANDS:")
# not sure how to identify these

print("\nThe rules for changes are found HERE")

# the run log does not currently feature the rule changes, only the initial terminal output does
# probably have to do it manually 

print("\nA list of all atoms that are modified can be found HERE")

# what's the best way to match residues from the original PDB file to their corresponding step1_out files? Coordinates?

print("\nWe have topology files for these ligands:") # how to check?

print("\nWe do not have topology files for these ligands:")

print("\nYou can (1) remove them from the input pdb file; (2) run with all atoms with zero charge; (3) make a topology file using instructions in:")

print("\n\n")
