#!/bin/bash

# Output file
output_file="bound_waters.txt"

# Write header to the output file
echo -e "Dir\tstep1_out\tfort.38" > "$output_file"

# Loop through each subdirectory
for dir in */; do
    # Ensure it's a directory
    if [ -d "$dir" ]; then
        # Initialize counts
        hoh_o_count=0
        hohdm_count=0
        
        # Count occurrences in step1_out.pdb
        if [ -f "$dir/step1_out.pdb" ]; then
            hoh_o_count=$(grep -Ec 'O   HOH' "$dir/step1_out.pdb")
        fi
        
        # Count occurrences in fort.38
        if [ -f "$dir/fort.38" ]; then
            hohdm_count=$(grep -c 'HOHDM' "$dir/fort.38")
        fi
        
        # Append results to output file
        echo -e "${dir%/}\t$hoh_o_count\t$hohdm_count" >> "$output_file"
    fi
done

echo "Results saved in $output_file"

