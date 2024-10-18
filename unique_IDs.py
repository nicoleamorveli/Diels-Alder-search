# Description:
# This script reads a text file containing IDs, identifies IDs that appear only once,
# and writes these unique IDs to an output text file.

# Inputs:
# - input_file_path: Path to the input text file containing IDs (one ID per line).
# - output_file_path: Path to the output text file where the unique IDs will be saved.

# Output:
# - A single text file containing IDs that appear only once in the input file, 
#   with each unique ID on a new line.

# Example:
# If you have the following input file (input.txt):
# A123
# B456
# C789
# A123
# D012

# Running the script as follows:
# python script.py input.txt output.txt

# The output.txt will contain:
# B456
# C789
# D012

import sys
from collections import Counter

def find_unique_ids(input_file_path, output_file_path):
    # Read the IDs from the text file
    with open(input_file_path, 'r') as file:
        ids = file.read().splitlines()

    # Count the occurrences of each ID
    id_counts = Counter(ids)

    # Filter IDs that appear only once
    unique_ids = [id for id, count in id_counts.items() if count == 1]

    # Write the unique IDs to the output file
    with open(output_file_path, 'w') as output_file:
        for id in unique_ids:
            output_file.write(f"{id}\n")

    print(f"Unique IDs have been saved to {output_file_path}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    find_unique_ids(input_file_path, output_file_path)

if __name__ == "__main__":
    main()
