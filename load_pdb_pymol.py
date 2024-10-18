# Description:
# This script generates PyMOL load commands for all PDB files located in a specified directory.
# It creates a script file that can be executed in PyMOL to load all PDB structures at once.

# Inputs:
# - A directory path containing PDB files: The script scans this directory for files ending with ".pdb".

# Output:
# - A script file named "load_pdb_files.pml" that contains load commands for each PDB file found in the directory.

# Example:
# If the directory contains the following PDB files:
# - protein1.pdb
# - protein2.pdb
# - protein3.pdb

# Running the script will produce a load_pdb_files.pml file containing:
# load /Users/nicolemorveli/Documents/Personal-project/SSN/small-networks/pdb_files/protein1.pdb
# load /Users/nicolemorveli/Documents/Personal-project/SSN/small-networks/pdb_files/protein2.pdb
# load /Users/nicolemorveli/Documents/Personal-project/SSN/small-networks/pdb_files/protein3.pdb

import os

# Path to the directory containing the PDB files
pdb_directory = "/Users/nicolemorveli/Documents/Personal-project/SSN/small-networks/pdb_files"

# Generate load commands for each PDB file
load_commands = []
for file_name in os.listdir(pdb_directory):
    if file_name.endswith(".pdb"):
        pdb_path = os.path.join(pdb_directory, file_name)
        load_commands.append(f"load {pdb_path}")

# Save the load commands to a script file
with open("load_pdb_files.pml", "w") as f:
    f.write("\n".join(load_commands))
