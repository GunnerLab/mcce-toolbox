import os
import pandas as pd
import re

def process_run_logs():
    # Get the base directory from the user
    base_dir = input("Enter the base directory: ").strip()

    # Ensure the base directory exists
    if not os.path.isdir(base_dir):
        print(f"The directory '{base_dir}' does not exist.")
        return

    # Create a "runs" directory within the base directory
    runs_dir = os.path.join(base_dir, "runs")
    os.makedirs(runs_dir, exist_ok=True)

    # Prepare a dataframe to store the results
    step_timings = pd.DataFrame()

    # Walk through all subdirectories
    for root, dirs, files in os.walk(base_dir):
        if "run.log" in files:
            run_log_path = os.path.join(root, "run.log")

            # Extract the parent directory name
            parent_dir = os.path.basename(root)

            # Initialize a dictionary to hold the timings for this run.log
            timing_data = {}

            # Read and process the run.log file
            with open(run_log_path, "r") as file:
                for line in file:
                    match_step1 = re.search(r"Total time for step1.*?([\d.]+)\s*seconds", line)
                    match_step2 = re.search(r"Total time for step2.*?([\d.]+)\s*seconds", line)
                    match_step3 = re.search(r"Total time for step3.*?([\d.]+)\s*seconds", line)
                    match_mc = re.search(r"Total time on MC:.*?([\d.]+)\s*seconds", line)

                    if match_step1:
                        timing_data["Step_1"] = float(match_step1.group(1))
                    elif match_step2:
                        timing_data["Step_2"] = float(match_step2.group(1))
                    elif match_step3:
                        timing_data["Step_3"] = float(match_step3.group(1))
                    elif match_mc:
                        timing_data["Step_4"] = float(match_mc.group(1))

            # Add the timings to the dataframe
            if timing_data:
                step_timings[parent_dir] = pd.Series(timing_data)

    # Transpose the dataframe so that steps are rows and directories are columns
    step_timings = step_timings.T
    step_timings.sort_index(ascending=True, inplace=True) # make sure proteins are in alphabetical order
    step_timings["Total"] = step_timings["Step_1"] + step_timings["Step_2"] + step_timings["Step_3"] + step_timings["Step_4"]

    # Export the dataframe as a .txt file
    output_file = os.path.join(runs_dir, f"{os.path.basename(base_dir)}_step_timings.txt")
    step_timings.to_csv(output_file, sep="\t")

    print(f"Step timings have been saved to {output_file}")

if __name__ == "__main__":
    process_run_logs()

