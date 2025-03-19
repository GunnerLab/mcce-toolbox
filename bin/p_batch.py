#!/usr/bin/env python
import os
import sys
import requests
import shutil
import json
from pathlib import Path

default_script_content = """#!/bin/bash
echo "Running default shell script in $(pwd)"

step1.py prot.pdb --dry -d 8
step2.py -d 8
step3.py -d 8
step4.py --xts
"""

# uses grafitti font from https://patorjk.com/software/taag/#p=display&f=Graffiti&t=Bench%20Batch%20v.1
ascii_art_open = r"""
__________                __         .__         __________         __         .__               ________
\______   \_______  _____/  |_  ____ |__| ____   \______   \_____ _/  |_  ____ |  |__   ___  __  \_____  \
 |     ___/\_  __ \/  _ \   __\/ __ \|  |/    \   |    |  _/\__  \\   __\/ ___\|  |  \  \  \/ /   /  ____/
 |    |     |  | \(  <_> )  | \  ___/|  |   |  \  |    |   \ / __ \|  | \  \___|   Y  \  \   /   /       \
 |____|     |__|   \____/|__|  \___  >__|___|  /  |______  /(____  /__|  \___  >___|  /   \_/ /\ \_______ \
                                   \/        \/          \/      \/          \/     \/        \/         \/
"""

def create_default_script():
    """Creates a shell script to perform a dry, level 1, d=8 run with delphi.
    Only run if a custom script is not given.
    """
    script_path = "default_script.sh"
    with open(script_path, "w") as script_file:
        script_file.write(default_script_content)
    os.chmod(script_path, 0o755)  # Make the script executable
    return script_path

def modify_script_for_runprm(script_path):
    """If run.prm.custom exists in PWD, modifies shell script to incorporate run.prm.custom.
    
    Args:
    script_path -- the shell script containing the base instructions for the MCCE run.
    """
    if os.path.exists("run.prm.custom"):
        with open(script_path, "r") as script_file:
            lines = script_file.readlines()
        with open(script_path, "w") as script_file:
            for line in lines:
                if any(cmd in line for cmd in ["step1.py", "step2.py", "step3.py", "step4.py"]) and " -load_runprm run.prm.custom" not in line:
                    line = line.strip() + " -load_runprm run.prm.custom\n"
                    script_file.write(line)
                else:
                    script_file.write(line) # hacky but will hopefully stop this function from wiping out the shell script
   
def fetch_protein(pdb_id):
    """Fetches a protein file from rcsb.org, creates a directory named after the protein, and symbolically links it to prot.pdb.
    
    Args:
    pdb_id -- four-character code corresponding to the protein in the database
    """
    pdb_id = pdb_id.lower()
    
    if len(pdb_id) == 4: # only make server request if pdb code is possible
        url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        response = requests.get(url)
    else:
        print("\nCode must be four letters!")
        return False

    if response.status_code == 200:
        os.makedirs(pdb_id, exist_ok=True)
        pdb_path = os.path.join(pdb_id, f'{pdb_id}.pdb')
        with open(pdb_path, 'wb') as file:
            file.write(response.content)
        
        link_path = os.path.join(pdb_id, 'prot.pdb')
        if os.path.exists(link_path) or os.path.islink(link_path):
            os.remove(link_path)
        os.symlink(pdb_path, link_path)
        
        print(f"\nProtein {pdb_id.upper()} downloaded successfully.")
        return True
    else:
        print(f"\nProtein {pdb_id.upper()} not found in RCSB database.")
        return False

def get_random_protein_codes(protein_code, num_samples=10):
    """Returns up to `num_samples` random similar protein codes for a given entry in the JSON file.

    Args: 
    protein_code -- what protein code are we getting similar structs for
    num_samples -- how many sim_structs do we want to compare against?
    """

    script_dir = os.path.dirname(os.path.abspath(__file__))
    json_filename = os.path.join(script_dir, "sim_structs.json")

    if not os.path.exists(json_filename):
        print("Error: JSON file not found.")
        return []

    with open(json_filename, 'r') as json_file:
        try:
            data = json.load(json_file)
        except json.JSONDecodeError:
            print("Error: JSON file is corrupted or empty.")
            return []

    if protein_code not in data:
        print(f"Protein code {protein_code} not found in the dataset.")
        return []

    similar_proteins = data[protein_code]
    return random.sample(similar_proteins, min(num_samples, len(similar_proteins)))

def process_protein_file(protein_path, script_path):
    """Executes MCCE scripts on dir containing a protein path.
    Makes sure all dirs contain their associated PDB file and run.prm.custom (if exists).

    Args:
    protein_path -- path to the protein to be run
    script_path -- path to the script 
    """

    protein_name = os.path.splitext(os.path.basename(protein_path))[0]
    protein_dir = os.path.abspath(protein_name)

    # Create directory for the protein
    os.makedirs(protein_dir, exist_ok=True)

    # Copy the protein file into the directory
    shutil.copy2(protein_path, os.path.join(protein_dir, os.path.basename(protein_path)))

    # Create a symbolic link to "prot.pdb"
    prot_pdb_path = os.path.join(protein_dir, "prot.pdb")
    if not os.path.exists(prot_pdb_path):
        os.symlink(os.path.basename(protein_path), prot_pdb_path)

    modify_script_for_runprm(script_path)

    # Create a symbolic link to "run.prm.custom" if it exists
    run_prm_path = os.path.join(os.getcwd(), "run.prm.custom")
    if os.path.exists(run_prm_path):
        run_prm_link = os.path.join(protein_dir, "run.prm.custom")
        if not os.path.exists(run_prm_link):
            os.symlink(run_prm_path, run_prm_link)

    # Execute the shell script in the directory in parallel without output
    os.system(f"cd {protein_dir} && bash ../{script_path} > /dev/null 2>&1 &")
   
def should_process_protein(protein_name):
    """Reads book.txt for c and x to check if should process this protein."""
    if os.path.exists("book.txt"):
        with open("book.txt", "r") as book_file:
            for line in book_file:
                if line.strip().startswith(protein_name):
                    if "c" in line or "x" in line:
                        return False
    return True

def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:

        print("Usage: p_batch.py <protein_files_or_directory> [<shell_script>]\n")

        print("Use flag '-h' or '-help' for more details.")
        sys.exit(1)

    if sys.argv[1] == "--help" or sys.argv[1] == "-help" or sys.argv[1] == "-h" or sys.argv[1] == "--h":

        print("\nProtein Batch accepts a file containing paths to protein files or a directory containing protein files. A shell script giving custom instructions for the bench batch may also be given. If custom instructions are not given, a default shell script will be created and executed.\n")

        print("Protein Batch creates a file called book.txt listing working proteins. If book.txt exists prior to execution, the program will read book.txt to know what runs to perform. In this way, a user may edit book.txt to run batches on a subset of desired proteins. File corresponding to lines containing ' c' or ' x' will not be run during a bench batch.\n")

        print("Usage: p_batch.py <protein_files_or_directory> [<shell_script>]")

        sys.exit(1)

    if sys.argv[1] == "rcsb":

        prot_code = input("Please input a four letter RCSB protein code (e.g. 4lzt): ").strip().lower()

        # maybe add a portion here if 4lzt already exists in a directory?

        if fetch_protein(prot_code):

            print("Download successful!") # obv needs to be added to, can't stop here
            
            get_random_protein_codes(prot_code)

            sys.exit(1) 

        else:

            print("\nCode not found at RCSB.org. Aborting.")
            sys.exit(1)

    else:

        input_path = sys.argv[1]
    
    script_path = sys.argv[2] if len(sys.argv) == 3 else create_default_script()
    
    # if file named custom.run.prm exists in the present working directory, make a symbolic link to it in each protein directory, and add "-load_rpm custom.run.prm" to each step. 

    if not os.path.exists(script_path):
        print(f"Error: Shell script '{script_path}' does not exist.")
        sys.exit(1)

    print(ascii_art_open)

    modify_script_for_runprm(script_path) # if run.prm.custom exists, make sure shell script sees it
    existing_dirs = []

    # Process individual protein files or all files in a directory
    if os.path.isfile(input_path): # parse individual protein files
        protein_name = os.path.splitext(os.path.basename(input_path))[0]
        print(protein_name)

        if os.path.exists(protein_name):
            existing_dirs.append(protein_name)

        process_protein_file(input_path, script_path)

    elif os.path.isdir(input_path): # parse a directory of protein files

        book_list = ""
        for filename in os.listdir(input_path):

            file_path = os.path.join(input_path, filename)

            if os.path.isfile(file_path) and filename.lower().endswith('.pdb'): # ONLY process PDB files
                protein_name = os.path.splitext(os.path.basename(file_path))[0]
                book_list += protein_name + "\n" # add to the book list

                if os.path.exists(protein_name):
                    existing_dirs.append(protein_name)
        
        if Path("book.txt").is_file(): # if we already have a book.txt

            print("\nbook.txt found! Protein files identified in book.txt: ")
            # print all valid files found in book.txt
            with open("book.txt") as f:
                reference_book = f.read() # only print lines w/o "c" and "x"
            print("\n" + reference_book)
            prots_to_run = ""
            for line in reference_book.split("\n"):
                if " c" not in line and " x" not in line:
                    prots_to_run += line + "\n"
                    print(line)

            if prots_to_run:
                print("These proteins will be run:\n\n" + prots_to_run +  "Pre-existing directories for these proteins will be emptied and replaced with information from the new run. ")

            else:
                print("No runnable proteins found. Clean the book and try again. Aborting MCCE.")
                sys.exit(1)

            if Path("run.prm.custom").is_file():
                print("\nrun.prm.custom found! The given shell script will be overwritten to read from run.prm.custom.")

            # LIST ALL SETTINGS, run.prm, extra.tpl, shell script, directories to be used, etc.
            response = input("Run MCCE with the current settings? (yes/y/no)").strip().lower() 
            
            # should empty protein file prior to processing so mcce_stat works properly

            if response == "yes" or response == "y":    
                
                modify_script_for_runprm(script_path)
                for filename in os.listdir(input_path):

                    if should_process_protein(filename[0:-4]) == True: # is filename minus ".pdb" in book.txt?
                        
                        file_path = os.path.join(input_path, filename)
                        if os.path.isfile(file_path) and filename.lower().endswith('.pdb'):
                            process_protein_file(file_path, script_path) # if so, process as usual
                            print("Processing " + file_path + "...")

                print("\nBash script is being executed in each directory. You can double check processes are being executed by running command 'top', or by running 'current_progress.py'")
            
            else:
                print("\nAborting MCCE")
                sys.exit(1)

        else: 
            # write to the book list here
            with open('book.txt', 'w') as file:
                file.write(book_list)

            print("\nNew book.txt created. You can remove protein files to be run by editing book.txt if desired, and resume by running bb again. ")

            print("These proteins will be run:\n\n" + book_list +  "Pre-existing directories for these proteins will be emptied and replaced with information from the new run. ")
            
            if Path("run.prm.custom").is_file():
                print("\nrun.prm.custom found! The given shell script will be overwritten to read from run.prm.custom.")

            response = input("Run MCCE with the current settings? (yes/y/no)").strip().lower()
            if response == "yes" or response == "y":
                modify_script_for_runprm(script_path)
                for filename in os.listdir(input_path): # check if filename is also in book.txt 
                    file_path = os.path.join(input_path, filename)
                    if os.path.isfile(file_path) and filename.lower().endswith('.pdb'):
                        process_protein_file(file_path, script_path)
                        print("Processing " + file_path + "...")

                print("\nBash script is being executed in each directory. You can double check processes are being executed by running command 'top', or by running 'current_progress.py'")
            else:
                print("\nAborting MCCE")
                sys.exit(1)

    else:

        print(f"\nError: '{input_path}' is neither a file nor a directory.")
        sys.exit(1)

if __name__ == "__main__":
    main()

