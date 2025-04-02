#!/usr/bin/env python
"""
Created on Mar 16 06:00:00 2025

@author: Gehan Ranepura
"""

import os
import sys

def remove_bk_lines(input_dir, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Loop through all files in the input directory
    for filename in os.listdir(input_dir):
        # Check if the file is a regular file
        if os.path.isfile(os.path.join(input_dir, filename)):
            input_file_path = os.path.join(input_dir, filename)
            output_file_path = os.path.join(output_dir, filename)

            try:
                with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
                    # Read each line in the input file
                    for line in infile:
                        # Write the line to the output file only if it doesn't contain 'BK'
                        if 'BK' not in line:
                            outfile.write(line)
                print(f"Processed {filename}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

def main():
    # Check if input and output directories are provided
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_directory> <output_directory>")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]

    # Call the function to remove BK lines
    remove_bk_lines(input_directory, output_directory)

if __name__ == "__main__":
    main()

