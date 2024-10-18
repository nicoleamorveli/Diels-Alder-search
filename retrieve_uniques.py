# Description:
# This script reads multiple input text files containing IDs and retrieves only the IDs 
# from the target file that do not appear in any of the comparison files. 
# It writes the unique IDs to an output text file.

# Inputs:
# - target_file: Path to the text file containing IDs to be compared.
# - comparison_files: List of paths to the text files for comparison (one or more).
# - output_file: Path to the output text file where the unique IDs will be written.

# Output:
# - A single text file containing the unique IDs from the target file, which do not appear in 
#   any of the comparison files, with each ID on a new line.

# Example:
# If you have the following files:
# target.txt:
# A123
# B456
# C789

# comparison1.txt:
# B456
# D012

# comparison2.txt:
# C789
# E345

# Running the script as follows:
# python script.py target.txt comparison1.txt comparison2.txt output.txt

# The output.txt will contain:
# A123

import sys

def retrieve_unique_ids(target_file, comparison_files, output_file):
    with open(target_file, 'r') as file:
        target_ids = set(line.strip() for line in file)

    comparison_ids = set()
    for file_path in comparison_files:
        with open(file_path, 'r') as file:
            comparison_ids.update(line.strip() for line in file)

    unique_ids = target_ids - comparison_ids

    with open(output_file, 'w') as output:
        for unique_id in sorted(unique_ids):
            output.write(f"{unique_id}\n")

def main():
    if len(sys.argv) < 4:
        print("Usage: python script.py <target_file> <comparison_file1> <comparison_file2> ... <output_file>")
        sys.exit(1)

    target_file = sys.argv[1]
    comparison_files = sys.argv[2:-1]
    output_file = sys.argv[-1]

    retrieve_unique_ids(target_file, comparison_files, output_file)

if __name__ == "__main__":
    main()
