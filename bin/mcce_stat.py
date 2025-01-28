#!/usr/bin/env python
import os

def check_files_in_directories():
    # Get the current working directory
    current_dir = os.getcwd()
    
    # List all subdirectories in the current directory, sorted alphabetically
    subdirectories = sorted([d for d in os.listdir(current_dir) if os.path.isdir(d)])

    # Define the target files
    target_files = ["step1_out.pdb", "step2_out.pdb", "head3.lst", "pK.out"]
    step_names = ["step1", "step2", "step3", "step4"]

    # Store results for output
    results = []
    all_pk_out_exists = True
    missing_step1_dirs = []

    for subdir in subdirectories:
        subdir_path = os.path.join(current_dir, subdir)

        # Check for the existence of the four target files
        files_status = {file: os.path.isfile(os.path.join(subdir_path, file)) for file in target_files}

        # Check timestamps if pK.out and head3.lst both exist
        if files_status.get("pK.out") and files_status.get("head3.lst"):
            pk_out_path = os.path.join(subdir_path, "pK.out")
            head3_path = os.path.join(subdir_path, "head3.lst")
            if os.path.getmtime(pk_out_path) < os.path.getmtime(head3_path):
                raise Exception(f"pK.out in {subdir} was created before head3.lst. head3.lst should always be created before pK.out. There may be a run in progress on a reused directory, or the run may be invalid and its result nonsense.")

        # Record the row for this directory
        results.append([subdir] + ["Exists" if files_status[file] else "Missing" for file in target_files])

        # Track directories missing step1_out.pdb
        if not files_status.get("step1_out.pdb"):
            missing_step1_dirs.append(subdir)

        # Check for conformer errors in head3.lst
        if files_status.get("head3.lst"):
            head3_path = os.path.join(subdir_path, "head3.lst")
            conf_error_lines = []

            try:
                with open(head3_path, 'r') as head3_file:
                    for line in head3_file:
                        columns = line.split()
                        if columns and columns[-1] == 'f':
                            conf_error_lines.append(line.strip())

                if conf_error_lines:
                    # Save the error lines to conf_errors.txt in the directory
                    error_file_path = os.path.join(subdir_path, "conf_errors.txt")
                    with open(error_file_path, 'w') as error_file:
                        error_file.write("\n".join(conf_error_lines))

            except Exception as e:
                print(f"Error reading head3.lst in {subdir}: {e}")
                return

        # Update all_pk_out_exists flag
        if not files_status.get("pK.out"):
            all_pk_out_exists = False

    # Print results in column format with proper spacing
    col_widths = [max(len(row[i]) for row in (["Directory"] + target_files) + results) for i in range(len(target_files) + 1)]
    header = ["Prot"] + step_names
    print(" ".join(header[i].ljust(col_widths[i]) for i in range(len(header))))
    for row in results:
        print(" ".join(row[i].ljust(col_widths[i]) for i in range(len(row))))

    # Report directories missing step1_out.pdb
    if missing_step1_dirs:
        print("These directories failed to start correctly: " + ", ".join(sorted(missing_step1_dirs)))

    # Final summary messages
    if all_pk_out_exists and len(results) == len(subdirectories):
        print("All runs completed as expected.")
    elif len(results) == len(subdirectories):
        print("All runs progressing as expected.")

if __name__ == "__main__":
    check_files_in_directories()

