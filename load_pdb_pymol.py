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
