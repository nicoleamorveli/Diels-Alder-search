import sys
import os

# Description:
# This script reads multiple input text files containing IDs (one per line),
# combines these IDs, removes duplicates, and writes the unique IDs to an output text file.
# It is intended to process lists of IDs from Foldseek search results.

# Inputs:
# - One or more input text files: Each file should contain a list of IDs (one ID per line).
# - The last argument is the output text file: This is where the unique IDs will be written.

# Output:
# - A single text file containing the unique IDs from the input files, with each ID on a new line.

# Example:
# If you have the following files:
# file1.txt:
# A123
# B456
# C789

# file2.txt:
# B456
# D012
# E345

# Running the script as follows:
# python combine_IDs.py file1.txt file2.txt output.txt

# The output.txt will contain:
# A123
# B456
# C789
# D012
# E345

def read_ids_from_files(input_files):
    ids = set()
    for input_file in input_files:
        if os.path.exists(input_file):
            with open(input_file, 'r') as file:
                ids.update(id.strip() for id in file.readlines())
        else:
            print(f"Warning: {input_file} does not exist and will be skipped.")
    return list(ids)

def write_unique_ids(output_file, unique_ids):
    with open(output_file, 'w') as file:
        for id in unique_ids:
            file.write(f"{id}\n")

def main(input_files, output_file):
    ids = read_ids_from_files(input_files)
    write_unique_ids(output_file, ids)
    print(f"Unique IDs have been written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python combine_IDs.py <input_file1.txt> <input_file2.txt> ... <output_file.txt>")
    else:
        input_files = sys.argv[1:-1]
        output_file = sys.argv[-1]
        main(input_files, output_file)

