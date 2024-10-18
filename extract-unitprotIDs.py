import sys
import csv

# Description:
# This script processes multiple Foldseek web result files in the m8 format, 
# extracts UniProt IDs from the second column of each file, and writes these IDs to an output CSV file.

# Inputs:
# - One or more input m8 files: Each file should be a result from Foldseek searches against the Alphafold proteome database, AFDB50, and SwissProt databases.
# - The last argument is the output CSV file: This is where the extracted UniProt IDs will be written.

# Output:
# - A single CSV file containing the extracted UniProt IDs from the input files, with each ID on a new line.

# Example:
# If you have the following files:
# file1.m8:
# query1   tr|A0A024R161-1|Uncharacterized protein OS=Homo sapiens OX=9606 GN=XXXX PE=1 SV=1
# query2   sp|P62258-2|TBB1_HUMAN Tubulin beta chain OS=Homo sapiens OX=9606 GN=TUBB PE=1 SV=2

# file2.m8:
# query3   sp|Q9Y6K9-3|H2B1M_HUMAN Histone H2B type 1-M OS=Homo sapiens OX=9606 GN=H2BC12 PE=1 SV=3
# query4   tr|A0A024R161-1|Uncharacterized protein OS=Homo sapiens OX=9606 GN=XXXX PE=1 SV=1

# Running the script as follows:
# python script.py file1.m8 file2.m8 output.csv

# The output.csv will contain:
# A0A024R161
# P62258
# Q9Y6K9
# A0A024R161


def extract_identifiers(input_files, output_file):
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        for input_file in input_files:
            with open(input_file, 'r') as infile:
                reader = csv.reader(infile, delimiter='\t')
                for line in reader:
                    if len(line) > 1:  # Ensure the line has at least two columns
                        identifier = line[1].split('-')[1]  # Extract the identifier
                        writer.writerow([identifier])

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py input1.m8 input2.m8 ... output.csv")
        sys.exit(1)

    input_files = sys.argv[1:-1]
    output_file = sys.argv[-1]

    extract_identifiers(input_files, output_file)
    print("Conversion completed.")

if __name__ == "__main__":
    main()

