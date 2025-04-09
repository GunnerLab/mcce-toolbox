#!/usr/bin/env python
import os
import csv
import re

def txt_to_csv(input_file):
    # converts pK.out into pK.csv for easy data handling
    output_file = input_file.rsplit('.', 1)[0] + ".csv"
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            stripped = line.lstrip()  # Remove leading spaces/tabs
            first_split = stripped.split(maxsplit=1)  # Split into two parts at the first whitespace

            if len(first_split) > 1:
                first_part, rest = first_split
                rest = ','.join(rest.split())  # Replace remaining spaces/tabs with commas
                outfile.write(f"{first_part},{rest}\n")
            else:
                outfile.write(stripped + '\n')

def collect_pk_files(root_dir):
    """Recursively create all 'pK.out.csv' files in subdirectories."""
    pk_files = []
    for subdir, _, files in os.walk(root_dir):
        if "pK.out" in files:
            txt_file = os.path.join(subdir, "pK.out")
            os.chdir(subdir)  # Change to the directory containing the file
            txt_to_csv(txt_file)  # Convert to CSV
            os.chdir(root_dir)  # Return to root directory
            pk_files.append(os.path.join(subdir, "pK.csv"))
    return pk_files

def merge_csv_files(pk_files, output_file):
    """Merge all 'pK.csv' files into one, keeping only the first file's header."""
    if not pk_files:
        print("No 'pK.csv' files found.")
        return

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    header_written = False
    rows = []
    
    with open(output_file, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv)
        
        for file in pk_files:
            directory = os.path.basename(os.path.dirname(file))  # Get the directory name
            
            with open(file, 'r', newline='') as in_csv:
                reader = csv.reader(in_csv)
                header = next(reader)  # Read header
                
                if not header_written:
                    writer.writerow(["Source Directory"] + header)  # Add source directory to header
                    rows.append(["Source Directory"] + header)
                    header_written = True
                
                header_length = len(["Source Directory"] + header)

                for row in reader:
                    full_row = [directory] + row
                    full_row += [""] * (header_length - len(full_row))  # Pad shorter rows
                    writer.writerow(full_row)
                    rows.append(full_row)

    txt_output_file = output_file.replace(".csv", ".txt")
    max_lengths = [max(len(str(entry)) for entry in column) for column in zip(*rows)]

    with open(txt_output_file, 'w') as txt_file:
        for row in rows:
            formatted_row = "  ".join(str(entry).ljust(max_lengths[i]) for i, entry in enumerate(row))
            txt_file.write(formatted_row + "\n")

    print(f"Formatted TXT file created at: {txt_output_file}")

if __name__ == "__main__":
    root_directory = os.getcwd()  # Change this to your target directory if needed
    output_directory = os.path.join(root_directory, "batch_analysis")
    output_file = os.path.join(output_directory, "all_pK_data.csv")
    
    pk_files = collect_pk_files(root_directory)
    merge_csv_files(pk_files, output_file)
    
    print(f"Merged CSV created at: {output_file}")
