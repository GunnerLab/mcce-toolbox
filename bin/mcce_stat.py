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
    step_weights = {"step1_out.pdb": 0.5, "step2_out.pdb": 0.5, "head3.lst": 7, "pK.out": 2}  # Step completion percentages

    # Store results for output
    results = []
    all_pk_out_exists = True
    missing_step1_dirs = []
    completed_count = 0
    completed_files = set()
    total_completion = 0

    for subdir in subdirectories:
        subdir_path = os.path.join(current_dir, subdir)

        # Only track directories containing "prot.pdb"
        if not os.path.isfile(os.path.join(subdir_path, "prot.pdb")):
            continue

        # Check for the existence of the four target files
        files_status = {file: os.path.isfile(os.path.join(subdir_path, file)) for file in target_files}

        # Check timestamps if pK.out and head3.lst both exist
        if files_status.get("pK.out") and files_status.get("head3.lst"):
            pk_out_path = os.path.join(subdir_path, "pK.out")
            head3_path = os.path.join(subdir_path, "head3.lst")
            if os.path.getmtime(pk_out_path) < os.path.getmtime(head3_path):
                raise Exception(f"pK.out in {subdir} was created before head3.lst")

        # Calculate step-wise completion
        protein_completion = sum(step_weights[file] for file in target_files if files_status[file])
        total_completion += protein_completion

        # Determine if step4 is completed (all files exist)
        step4_completed = all(files_status[file] for file in target_files)
        if step4_completed:
            completed_count += 1
            completed_files.add(subdir)

        # Record the row for this directory
        results.append([subdir] + ["Exists" if files_status[file] else "      " for file in target_files] + (["Done"] if step4_completed else [""]))

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

    # Write to book.txt
    book_path = os.path.join(current_dir, "book.txt")
    with open(book_path, 'w') as book_file:
        for subdir in subdirectories:
            if os.path.isfile(os.path.join(current_dir, subdir, "prot.pdb")):
                book_file.write(subdir + ("   c" if subdir in completed_files else "") + "\n")

    # Calculate completion percentage
    total_dirs = len(results)
    max_completion = total_dirs * sum(step_weights.values()) if total_dirs > 0 else 1  # Avoid division by zero
    completion_percentage = (total_completion / max_completion) * 100

    # Print percentage of completed rows
    print(f"Completion: {completion_percentage:.2f}%")

    # Print results in column format with proper spacing
    if results:
        col_widths = [max(len(str(row[i])) for row in ([["Directory"] + target_files + ["Status"]] + results)) for i in range(len(results[0]))]
    else:
        col_widths = [len(col) for col in ["Directory"] + target_files + ["Status"]]

    header = ["Directory"] + step_names + ["Status"]
    print(" ".join(header[i].ljust(col_widths[i]) for i in range(len(header))))
    for row in results:
        print(" ".join(str(row[i]).ljust(col_widths[i]) for i in range(len(row))))

    # Report directories missing step1_out.pdb
    if missing_step1_dirs:
        print("These directories failed to start correctly: " + ", ".join(sorted(missing_step1_dirs)))

    # Final summary messages
    if all_pk_out_exists and total_dirs == len(subdirectories):
        print("All runs completed as expected.")
    elif total_dirs == len(subdirectories):
        print("All runs progressing as expected.")

if __name__ == "__main__":
    check_files_in_directories()

