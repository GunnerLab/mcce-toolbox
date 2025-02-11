#!/usr/bin/env python
import os
import sys
import shutil
from pathlib import Path

default_script_content = """#!/bin/bash
echo "Running default shell script in $(pwd)"

step1.py prot.pdb --dry -d 8
step2.py -d 8
step3.py -d 8
step4.py --xts
"""

# uses grafitti font from https://patorjk.com/software/taag/#p=display&f=Graffiti&t=Bench%20Batch%20v.1
ascii_art_open = """
__________                     .__      __________         __         .__               ____ 
\______   \ ____   ____   ____ |  |__   \______   \_____ _/  |_  ____ |  |__   ___  __ /_   |
 |    |  _// __ \ /    \_/ ___\|  |  \   |    |  _/\__  \|_  __\/ ___\|  |  \  \  \/ /  |   |
 |    |   \  ___/|   |  \  \___|   Y  \  |    |   \ / __ \|  | \  \___|   Y  \  \   /   |   |
 |______  /\___  >___|  /\___  >___|  /  |______  /(____  /__|  \___  >___|  /   \_/ /\ |___|
        \/     \/     \/     \/     \/          \/      \/          \/     \/        \/      
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
    
def cleanup(d):
    """Used to delete occupied directories prior to running, to avoid confusion in case of errors and left-over data.
    
    Args:
    d -- a protein file, from which we remove the ".pdb" portion to access its directory.
    """

    trash_file = Path(d).stem
    shutil.rmtree(trash_file)

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

        print("Usage: bb.py <protein_files_or_directory> [<shell_script>]")
        sys.exit(1)

    if sys.argv[1] == "--help":

        print("\nBench Batch, or bb, accepts a file containing paths to protein files or a directory containing protein files. A shell script giving custom instructions for the bench batch may also be given. If custom instructions are not given, a default shell script will be created and executed.\n")

        print("bb creates a file called book.txt listing working proteins. If book.txt exists prior to bb being executed, bb will read book.txt to know what runs to perform. In this way, a user may edit book.txt to run batches on a subset of desired proteins.\n")

        print("Usage: bb.py <protein_files_or_directory> [<shell_script>]")

        sys.exit(1)

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

            if prots_to_run:
                print("These proteins will be run:\n\n" + prots_to_run +  "Pre-existing directories for these proteins will be emptied and replaced with information from the new run. ")

            else:
                print("No runnable proteins found. Clean the book and try again. Aborting MCCE.")
                sys.exit(1)

            if Path("run.prm.custom").is_file():
                print("\nrun.prm.custom found! The given shell script will be overwritten to read from run.prm.custom.")

            # LIST ALL SETTINGS, run.prm, extra.tpl, shell script, directories to be used, etc.
            response = input("Run MCCE with the current settings? (yes/y)").strip().lower() 
            
            # should empty protein file prior to processing so mcce_stat works properly

            if response == "yes" or response == "y":    
                
                modify_script_for_runprm(script_path)
                for filename in os.listdir(input_path):

                    if should_process_protein(filename[0:-4]) == True: # is filename minus ".pdb" in book.txt?
                        
                        file_path = os.path.join(input_path, filename)
                        if os.path.isfile(file_path):
                            # cleanup(file_path)
                            process_protein_file(file_path, script_path) # if so, process as usual
                            print("Processing " + file_path + "...")

        else: 
            # write to the book list here
            with open('book.txt', 'w') as file:
                file.write(book_list)

            print("\nNew book.txt created. You can remove protein files to be run by editing book.txt if desired, and resume by running bb again. ")

            print("These proteins will be run:\n\n" + book_list +  "Pre-existing directories for these proteins will be emptied and replaced with information from the new run. ")
            
            if Path("run.prm.custom").is_file():
                print("\nrun.prm.custom found! The given shell script will be overwritten to read from run.prm.custom.")

            response = input("Run MCCE with the current settings? (yes/y)").strip().lower()
            if response == "yes" or response == "y":
                modify_script_for_runprm(script_path)
                for filename in os.listdir(input_path): # check if filename is also in book.txt 
                    file_path = os.path.join(input_path, filename)
                    if os.path.isfile(file_path):
                        # cleanup(filename) # presents issues at present
                        process_protein_file(file_path, script_path)
                        print("Processing " + file_path + "...")
            else:
                print("\nAborting MCCE")
                sys.exit(1)

    else:

        print(f"\nError: '{input_path}' is neither a file nor a directory.")
        sys.exit(1)

    print("\nBash script is being executed in each directory. You can double check processes are being executed by running command 'top', or by running 'current_progress.py'")

if __name__ == "__main__":
    main()

