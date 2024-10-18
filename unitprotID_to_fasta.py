import argparse
from Bio import SeqIO
from Bio import ExPASy

def retrieve_fasta(uniprot_ids, output_file):
    with open(output_file, "w") as f_out:
        for uniprot_id in uniprot_ids:
            try:
                with ExPASy.get_sprot_raw(uniprot_id) as handle:
                    record = SeqIO.read(handle, "swiss")
                    f_out.write(record.format("fasta"))
            except Exception as e:
                print(f"Error retrieving sequence for UniProt ID {uniprot_id}: {e}")

def main(input_file, output_file):
    with open(input_file, 'r') as f_in:
        uniprot_ids = [line.strip() for line in f_in if line.strip()]
    
    retrieve_fasta(uniprot_ids, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve FASTA sequences from UniProt IDs")
    parser.add_argument("input_file", help="Input text file containing UniProt IDs")
    parser.add_argument("output_file", help="Output text file to save FASTA sequences")
    args = parser.parse_args()

    main(args.input_file, args.output_file)

