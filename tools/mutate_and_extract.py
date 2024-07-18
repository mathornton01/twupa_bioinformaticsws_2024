import argparse
from Bio import SeqIO

def extract_and_modify_sequence(fasta_file, chromosome, location, nucleotide, bases_before, bases_after, output_file):
    # Read all sequences from the fasta file
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

    if chromosome not in sequences:
        raise ValueError(f"Chromosome {chromosome} not found in the provided FASTA file.")

    sequence = sequences[chromosome]

    # Adjust for 0-based indexing in Python
    location -= 1

    # Extract the required sequence
    start = max(0, location - bases_before)
    end = min(len(sequence), location + bases_after + 1)
    extracted_sequence = sequence[start:end]

    # Replace the nucleotide at the specified location
    if 0 <= location < len(sequence):
        modified_sequence = (
            extracted_sequence[:location - start] +
            nucleotide +
            extracted_sequence[location - start + 1:]
        )
    else:
        raise ValueError("Location is out of the range of the sequence")

    # Write the modified sequence to the output fasta file
    with open(output_file, "w") as out_file:
        out_file.write(f">{chromosome}_modified\n")
        for i in range(0, len(modified_sequence), 60):
            out_file.write(modified_sequence[i:i+60] + "\n")

def main():
    parser = argparse.ArgumentParser(description="Extract and modify a sequence from a FASTA file.")
    parser.add_argument("fasta_file", type=str, help="Path to the reference FASTA file.")
    parser.add_argument("chromosome", type=str, help="Chromosome identifier to extract the sequence from.")
    parser.add_argument("location", type=int, help="1-based index location of the nucleotide to be modified.")
    parser.add_argument("nucleotide", type=str, help="Nucleotide letter to replace at the given location.")
    parser.add_argument("bases_before", type=int, help="Number of bases to extract before the given location.")
    parser.add_argument("bases_after", type=int, help="Number of bases to extract after the given location.")
    parser.add_argument("output_file", type=str, help="Path to the output FASTA file.")

    args = parser.parse_args()

    extract_and_modify_sequence(
        args.fasta_file, 
        args.chromosome,
        args.location, 
        args.nucleotide, 
        args.bases_before, 
        args.bases_after, 
        args.output_file
    )

if __name__ == "__main__":
    main()

