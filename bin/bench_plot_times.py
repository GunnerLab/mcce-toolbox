#!/usr/bin/env python
"""
Created on Thu Oct 03 16:01:00 2024

@author: Gehan Ranepura
"""

import glob
import pandas as pd
import matplotlib.pyplot as plt

# Set data parameters for plot naming
pb_solver = "Delphi"
eps_in = 4
charge_set = "Parse"

# Collect runtime lines from run.log files in the the runs folder
log_files = glob.glob("runs/*/run.log")
search_terms = ["Total time for step1", "Total time for step2", "Total time for step3", "Total time on MC"]
with open("runtime_lines.txt", "w") as output_file:
    for log in log_files:
        with open(log, "r") as f:
            for line in f:
                if any(term in line for term in search_terms):
                    output_file.write(f"{log}: {line}")

# Read the file and initialize a dictionary to store the times
times = {
    'PDB Code': [],
    'Step 1 (premcce)': [],
    'Step 2 (rotamer making)': [],
    'Step 3 (energies)': [],
    'Step 4 (MC sampling)': []
}

# Open the file and parse each line
with open('runtime_lines.txt', 'r') as file:
    for line in file:
        if 'Total time for step1 (premcce)' in line:
            parts = line.split()
            pdb_code = parts[0].split('/')[-2]  # Extract PDB code
            time = float(parts[-2])/3600  # Extract time in hours
            times['PDB Code'].append(pdb_code)
            times['Step 1 (premcce)'].append(time)
            times['Step 2 (rotamer making)'].append(0)  # Initialize with 0
            times['Step 3 (energies)'].append(0)  # Initialize with 0
            times['Step 4 (MC sampling)'].append(0)  # Initialize with 0

        elif 'Total time for step2 (rotamer making)' in line:
            parts = line.split()
            pdb_code = parts[0].split('/')[-2]
            time = float(parts[-2])/3600
            if pdb_code in times['PDB Code']:
                idx = times['PDB Code'].index(pdb_code)
                times['Step 2 (rotamer making)'][idx] = time

        elif 'Total time for step3' in line:
            parts = line.split()
            pdb_code = parts[0].split('/')[-2]
            time = float(parts[-2])/3600
            if pdb_code in times['PDB Code']:
                idx = times['PDB Code'].index(pdb_code)
                times['Step 3 (energies)'][idx] = time

        elif 'Total time on MC' in line:
            parts = line.split()
            pdb_code = parts[0].split('/')[-2]
            time = float(parts[-2])/3600
            if pdb_code in times['PDB Code']:
                idx = times['PDB Code'].index(pdb_code)
                times['Step 4 (MC sampling)'][idx] = time

# Create a DataFrame from the times dictionary
df = pd.DataFrame(times)

# Set PDB Code as the index
df.set_index('PDB Code', inplace=True)

# Calculate the total time for sorting
df['Total Time'] = df.sum(axis=1)

# Sort the DataFrame by total time
df = df.sort_values(by='Total Time')

# Drop the 'Total Time' column for plotting
df = df.drop(columns='Total Time')

# Plotting the stacked bar graph
ax = df.plot(kind='bar', stacked=True, figsize=(12, 8))

# Adding the total time on top of each bar with 45-degree rotation
for i in range(len(df)):
    total_time = df.iloc[i].sum()
    ax.text(i, total_time, f'{total_time:.2f}', ha='center', va='bottom', rotation=45)

# Customizing the plot
plt.title(f'{pb_solver} Run Times for {charge_set} Level 2 Benchmarks (eps={eps_in})')
plt.xlabel('PDB Code')
plt.ylabel('Time (hours)')
plt.legend(title='Steps')
plt.xticks(rotation=45)
plt.tight_layout()

# Save the plot as a PNG file
plt.savefig(f'runtimePlot_{pb_solver}-d{eps_in}_{charge_set}.png', format='png', dpi=300)

# Display the plot
plt.show()

