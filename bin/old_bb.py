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

def create_default_script():
    script_path = "default_script.sh"
    with open(script_path, "w") as script_file:
        script_file.write(default_script_content)
    os.chmod(script_path, 0o755)  # Make the script executable
    return script_path

def process_protein_file(protein_path, script_path):
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

    # Execute the shell script in the directory in parallel without output
    os.system(f"cd {protein_dir} && bash ../{script_path} > /dev/null 2>&1 &")

def prompt_and_cleanup(existing_dirs):
    print("The following directories already exist:")
    for d in existing_dirs:
        print(f"- {d}")
    response = input("Do you want to re-run tests in the current directories? This will delete the current directories. (yes/no): ").strip().lower()
    if response == "yes" or response == "y":
        for d in existing_dirs:
            shutil.rmtree(d)
        return True
    else:
        print("Aborting.")
        sys.exit(1)

def clean_book(file_path): # used to clean book.txt, if desired
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        with open(file_path, 'w') as file:
            for line in lines:
                if line.strip().endswith('c'):
                    file.write(line.rstrip()[:-1] + '\n')
                else:
                    file.write(line)
        
        print("File cleaned successfully.")
    except FileNotFoundError:
        print(f"Error: {file_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:

        print("Usage: bb_v3.py <protein_files_or_directory> [<shell_script>]")
        sys.exit(1)

    if sys.argv[1] == "--help":

        print("\nBench Batch, or bb, accepts a file containing paths to protein files or a directory containing protein files. A shell script giving custom instructions for the bench batch may also be given. If custom instructions are not given, a default shell script will be created and executed.\n")

        print("bb creates a file called book.txt listing working proteins. If book.txt exists prior to bb being executed, bb will read book.txt to know what runs to perform. In this way, a user may edit book.txt to run batches on a subset of desired proteins.\n")

        print("Usage: bb_v3.py <protein_files_or_directory> [<shell_script>]")

        sys.exit(1)

    input_path = sys.argv[1]
    script_path = sys.argv[2] if len(sys.argv) == 3 else create_default_script()
    
    # if file named custom.run.prm exists in the present working directory, make a symbolic link to it in each protein directory, and add "-load_rpm custom.run.prm" to each step. 

    if not os.path.exists(script_path):
        print(f"Error: Shell script '{script_path}' does not exist.")
        sys.exit(1)

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

            if os.path.isfile(file_path):
                protein_name = os.path.splitext(os.path.basename(file_path))[0]
                book_list += protein_name + "\n" # add to the book list

                if os.path.exists(protein_name):
                    existing_dirs.append(protein_name)
        
        if Path("book.txt").is_file(): # if we already have a book.txt

            response = input("\nbook.txt found! Run protein files identified in book.txt? (yes/no): ").strip().lower()
            if response == "yes" or response == "y":
                with open("book.txt") as f:
                    reference_book = f.read().replace('\n', ' ')
                if existing_dirs:
                    if not prompt_and_cleanup(existing_dirs):
                        sys.exit(1)

                for filename in os.listdir(input_path):

                    if filename[0:-4] in reference_book: # is filename minus ".pdb" in book.txt?
                        
                        file_path = os.path.join(input_path, filename)
                        if os.path.isfile(file_path):
                            process_protein_file(file_path, script_path) # if so, process as usual
            
                print("\nBash script is being executed on directories existing in book.txt. You can double check processes are being executed by running command 'top', or by running 'current_progress.py'\n")

                sys.exit(1) # conclude program upon beginning processing

            else: 
                response = input("\nRun bb directly from given protein list/directory, disregarding current book.txt? (yes/no)")

                if response == "yes" or response == "y": 
                        
                    print("\nbook.txt disregarded, resume bb as normal.")

                else:

                    print("\nAborting.")

                    sys.exit(1)

            # should have some more error messages on here
                    
                if existing_dirs:
                    if not prompt_and_cleanup(existing_dirs):
                        sys.exit(1)
                for filename in os.listdir(input_path):
                    file_path = os.path.join(input_path, filename)
                    if os.path.isfile(file_path):
                        process_protein_file(file_path, script_path)
                        print("Processing " + file_path + "...")

        else: # get rid of else so if book.txt is disregarded we end up down here anyways
            # write to the book list here
            with open('book.txt', 'w') as file:
                file.write(book_list)

            response = input("\nNew book.txt created. You can remove protein files to be run by editing book.txt if desired, and resume by running bb again. Run all protein files immediately? (yes/no): ").strip().lower()
            if response == "yes" or response == "y":
                if existing_dirs:
                    if not prompt_and_cleanup(existing_dirs):
                        sys.exit(1)
                for filename in os.listdir(input_path): # check if filename is also in book.txt 
                    file_path = os.path.join(input_path, filename)
                    if os.path.isfile(file_path):
                        process_protein_file(file_path, script_path)
                        print("Processing " + file_path + "...")
            else:
                print("\nAborting. Edit book.txt to remove undesired proteins, and run bb again.")
                sys.exit(1)
        
    else:

        print(f"\nError: '{input_path}' is neither a file nor a directory.")
        sys.exit(1)

    print("\nBash script is being executed in each directory. You can double check processes are being executed by running command 'top', or by running 'current_progress.py'")

if __name__ == "__main__":
    main()

