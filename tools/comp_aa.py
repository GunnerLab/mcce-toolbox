import os
import pandas as pd

def merge_pka_files(dir1, dir2, output_file='merged_pkas.txt'):
    # Define paths to analysis folders
    analysis_path1 = os.path.join(dir1, "analysis")
    analysis_path2 = os.path.join(dir2, "analysis")

    # Check if analysis folders exist
    #if not os.path.isdir(analysis_path1):
    #    raise FileNotFoundError(f"The analysis folder is missing in the directory: {dir1}")
    #if not os.path.isdir(analysis_path2):
    #    raise FileNotFoundError(f"The analysis folder is missing in the directory: {dir2}")

    # Define paths to matched_pkas.txt files
    file_path1 = os.path.join(analysis_path1, "matched_pkas.txt")
    file_path2 = os.path.join(analysis_path2, "matched_pkas.txt")

    cwd = os.getcwd()
    print(cwd)

    # Check if matched_pkas.txt files exist
    #if not os.path.isfile(file_path1):
    #    raise FileNotFoundError(f"The file 'matched_pkas.txt' is missing in: {analysis_path1}")
    #if not os.path.isfile(file_path2):
    #    raise FileNotFoundError(f"The file 'matched_pkas.txt' is missing in: {analysis_path2}")
    
    # Load files as data frames
    df1 = pd.read_csv(file_path1, delim_whitespace=True)
    df2 = pd.read_csv(file_path2, delim_whitespace=True)

    df1.drop(columns=['expl'], inplace=True)
    
    # Drop the 'delta' column from each DataFrame
   # df1 = df1.drop(columns=['delta'], errors='ignore')
   # df2 = df2.drop(columns=['delta'], errors='ignore')

    # Merge the data frames on 'resid'
    merged_df = pd.merge(df1, df2, on='resid', suffixes=('_dir1', '_dir2'))

    column_to_move = merged_df.pop("expl")
    merged_df.insert(2, "expl", column_to_move)
    merged_df["change_between_directories"] = merged_df.iloc[:,6] - merged_df.iloc[:,3] 

    # Save the merged data frame to a .txt file
    merged_df.to_string(output_file, index=False)

    print(f"Merged file saved as: {output_file}")

print("This program asks for directories containing a folder named 'analysis', that contains a file named 'matched_pkas.txt'.\n")

# Example usage:
dir1 = input("Enter the path for the first directory: ")
dir2 = input("Enter the path for the second directory: ")
merge_pka_files(dir1, dir2)
