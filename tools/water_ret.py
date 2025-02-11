import os
import subprocess

def count_matches(pattern, file_path):
    """Use grep to count the number of lines matching the given pattern in the specified file."""
    try:
        result = subprocess.run(
                ['grep', f'^{pattern}', file_path, '|', 'wc', '-l'],
            # grep "^HOHDM" fort.38 | awk '{count++; sum+=$8} END {print "Rows:", count, "\nSum of 7.0 column:", sum, "\nWater Retention:", sum / count}'

            capture_output=True, text=True, shell=True, check=True
        )
        return int(result.stdout.strip())
    except subprocess.CalledProcessError:
        return 0

def main():
    print("Rate of Water Retention")
    for root, dirs, files in os.walk("."):
        if "fort.38" in files:
            file_path = os.path.join(root, "fort.38")
            
            # Get counts for waters and dm_waters
            waters = count_matches("HOH", file_path)
            dm_waters = count_matches("HOHDM", file_path)
            
            # Calculate the rate of water retention
            if waters > 0:
                retention_rate = (waters - dm_waters) / waters
            else:
                retention_rate = 0  # Avoid division by zero

            # Print the results for this directory
            print(f"Directory: {root}")
            print(f"waters: {waters}, dm_waters: {dm_waters}, retention_rate: {retention_rate:.4f}")
            print("-" * 40)

if __name__ == "__main__":
    main()

