#!/usr/bin/env python
import os
import sys
import shutil

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

    # Execute the shell script in the directory (hopefully in parallel)
    os.system(f"cd {protein_dir} && bash ../{script_path} &")

def prompt_and_cleanup(existing_dirs):
    print("The following directories already exist:")
    for d in existing_dirs:
        print(f"- {d}")
    response = input("Do you want to re-run tests in the current directories? This will delete the current directories. (yes/no): ").strip().lower()
    if response == "yes":
        for d in existing_dirs:
            shutil.rmtree(d)
        return True
    else:
        print("Aborting.")
        sys.exit(1)

def main():
    if len(sys.argv) != 3:
        print("Usage: python bench_bench_j.py <protein_files_or_directory> <shell_script>")
        sys.exit(1)

    input_path = sys.argv[1]
    script_path = sys.argv[2]

    if not os.path.exists(script_path):
        print(f"Error: Shell script '{script_path}' does not exist.")
        sys.exit(1)

    existing_dirs = []

    # Process individual protein files or all files in a directory
    if os.path.isfile(input_path):
        protein_name = os.path.splitext(os.path.basename(input_path))[0]
        if os.path.exists(protein_name):
            existing_dirs.append(protein_name)
        process_protein_file(input_path, script_path)
    elif os.path.isdir(input_path):
        for filename in os.listdir(input_path):
            file_path = os.path.join(input_path, filename)
            if os.path.isfile(file_path):
                protein_name = os.path.splitext(os.path.basename(file_path))[0]
                if os.path.exists(protein_name):
                    existing_dirs.append(protein_name)
        if existing_dirs:
            if not prompt_and_cleanup(existing_dirs):
                sys.exit(1)
        for filename in os.listdir(input_path):
            file_path = os.path.join(input_path, filename)
            if os.path.isfile(file_path):
                process_protein_file(file_path, script_path)
    else:
        print(f"Error: '{input_path}' is neither a file nor a directory.")
        sys.exit(1)

if __name__ == "__main__":
    main()

