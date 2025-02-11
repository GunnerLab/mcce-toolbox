import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def clean_pKa_column(column):
    """Cleans and converts pKa/Em column to float."""
    return column.str.replace(r"[><]", "", regex=True).astype(float)

def process_directories(dir1, dir2):
    # Locate the "pK.out" file in each directory
    file1 = os.path.join(dir1, "pK.out")
    file2 = os.path.join(dir2, "pK.out")
    
    # Ensure the files exist
    if not os.path.isfile(file1):
        raise FileNotFoundError(f"File 'pK.out' not found in {dir1}")
    if not os.path.isfile(file2):
        raise FileNotFoundError(f"File 'pK.out' not found in {dir2}")
    
    # Read the files into pandas DataFrames
    df1 = pd.read_csv(file1, delim_whitespace=True, usecols=["pH", "pKa/Em"])
    df2 = pd.read_csv(file2, delim_whitespace=True, usecols=["pH", "pKa/Em"])
    
    # Clean and convert the "pKa/Em" columns to float
    df1["pKa/Em"] = clean_pKa_column(df1["pKa/Em"].astype(str))
    df2["pKa/Em"] = clean_pKa_column(df2["pKa/Em"].astype(str))
    
    # Rename columns for clarity
    df1.rename(columns={"pKa/Em": "pKa/Em_dir1"}, inplace=True)
    df2.rename(columns={"pKa/Em": "pKa/Em_dir2"}, inplace=True)
    
    # Merge the two DataFrames on the "pH" column
    combined_df = pd.merge(df1, df2, on="pH")
    
    # Extract labels from the first three characters of "pH"
    combined_df['Label'] = combined_df['pH'].astype(str).str[:3]
    
    # Assign colors for unique labels
    unique_labels = combined_df['Label'].unique()
    colors = cm.tab10(range(len(unique_labels)))
    color_map = {label: colors[i] for i, label in enumerate(unique_labels)}
    
    # Plot the data
    plt.figure(figsize=(10, 8))
    for i, row in combined_df.iterrows():
        x, y = row["pKa/Em_dir1"], row["pKa/Em_dir2"]
        label = row['Label']
        plt.scatter(x, y, color=color_map[label], label=label if label not in plt.gca().get_legend_handles_labels()[1] else "", s=50)
        #plt.text(x, y, label, color=color_map[label], fontsize=9, ha='right', va='bottom')
    
    # Add reference line
    plt.plot([combined_df["pKa/Em_dir1"].min(), combined_df["pKa/Em_dir1"].max()], 
             [combined_df["pKa/Em_dir1"].min(), combined_df["pKa/Em_dir1"].max()], 
             color='black', linestyle='--', label='y=x')
   
    plt.plot([combined_df["pKa/Em_dir1"].min(), combined_df["pKa/Em_dir1"].max()],
             [combined_df["pKa/Em_dir1"].min() + 1, combined_df["pKa/Em_dir1"].max() + 1],
             color='grey', linestyle='--', label='y=x-1')

    plt.plot([combined_df["pKa/Em_dir1"].min(), combined_df["pKa/Em_dir1"].max()],
             [combined_df["pKa/Em_dir1"].min() - 1, combined_df["pKa/Em_dir1"].max() - 1],
             color='grey', linestyle='--', label='y=x+1')

    plt.plot([combined_df["pKa/Em_dir1"].min(), combined_df["pKa/Em_dir1"].max()],
             [combined_df["pKa/Em_dir1"].min() + 2, combined_df["pKa/Em_dir1"].max() + 2],
             color='grey', linestyle='--', label='y=x-2')

    plt.plot([combined_df["pKa/Em_dir1"].min(), combined_df["pKa/Em_dir1"].max()],
             [combined_df["pKa/Em_dir1"].min() - 2, combined_df["pKa/Em_dir1"].max() - 2],
             color='grey', linestyle='--', label='y=x+2')

    plt.xlabel(f"pKa/Em from {dir1}")
    plt.ylabel(f"pKa/Em from {dir2}")
    plt.title(f"Direct Comparison of pKa/Em Values, {dir1} v. {dir2}")
    plt.legend(title="Labels")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"output_pkas_{dir1}_v_{dir2}.png")

# Example usage:
if __name__ == "__main__":
    dir1 = input("Enter the path to the first directory: ")
    dir2 = input("Enter the path to the second directory: ")
    try:
        process_directories(dir1, dir2)
    except Exception as e:
        print(f"Error: {e}")
