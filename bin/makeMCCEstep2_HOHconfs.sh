#!/bin/bash
# Created on Dec 12 23:25:00 2024
# @author: Gehan Ranepura

# Display help message
function show_help() {
    echo "Usage: $0 <input_pdb> [-N <value>] [-h]"
    echo
    echo "Options:"
    echo "  -N <value>   Set the value for N (default is 25)"
    echo "  -h           Show this help message"
    echo
    echo "Description:"
    echo "  This script processes a PDB file, extracts water molecules, generates conformers, and"
    echo "  outputs a concatenated PDB file with the results."
    echo
    echo "Examples:"
    echo "  $0 input.pdb -N 5"
    echo "  $0 input.pdb"
    exit 0
}

# Check if the file argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input_pdb_file> [-N <value>] [-h]"
    exit 1
fi

# Set a variable called "input_pdb" to the first argument passed, Set default N value
input_pdb="$1"
N=25

# Shift the positional argument to process flags
shift

# Process the command line arguments for the -N and -h options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -N)
            # Ensure -N is followed by a valid number
            if [[ ! "$2" =~ ^[0-9]+$ ]]; then
                echo "Error: -N requires a numeric value."
                exit 1
            fi
            N="$2"
            shift 2
            ;;
        N)
            echo "Error: -N is correct option"
            exit 1
            ;; 
        -h|--help)
            show_help
            ;;
        *)
            echo "Error: Unrecognized argument '$1'."
            exit 1
            ;;
    esac
done

# Generate directory names and store them in a variable
dirs=$(awk '($4 != "" && $5 != "" && $6 != "") {print substr($4,1,3) "-" substr($5,1,1) "-" substr($6,1,3)}' "$input_pdb")
num_dirs=$(echo "$dirs" | wc -l)
echo "$N conformers will be made for the $num_dirs waters in $input_pdb"

# Remove any existing directories and recreate them
rm -rf $dirs
mkdir -p $dirs

# Create an array to store the paths of the newly created PDB files
new_pdb_files=()

# Process each newly created directory
for dir in $dirs; do
    echo "Processing water: $dir"

    # Ensure the directory exists before proceeding
    if [ ! -d "$dir" ]; then
        echo "Directory $dir does not exist. Skipping..."
        continue
    fi

    # Extract the water molecules corresponding to the directory name and write to 'wat.pdb'
    grep "$(echo "$dir" | sed 's,\-, ,g' | sed 's,\/,,g')" "$input_pdb" > "$dir/wat.pdb"

    # Run the Python script and write the output to 'HOH_coor.txt'
    pdb2mcce_step2HOHconfs.py "$dir/wat.pdb" -N "$N" > "$dir/HOH_coor.txt"

    # Move the generated PDB file to the current directory
    mv HOH_confs.pdb "$dir/"

    # Store the path to the generated PDB file for concatenation
    new_pdb_files+=("$dir/HOH_confs.pdb")
done

# Concatenate only the newly created PDB files
output_file="${input_pdb%.*}_step2_out.${input_pdb##*.}"
cat "${new_pdb_files[@]}" > "$output_file"

# Clean up the created directories
rm -rf $dirs

# Final success message
echo -e "Water conformers generated successfully, Congrats!!! \n The file '$output_file' has been generated."


