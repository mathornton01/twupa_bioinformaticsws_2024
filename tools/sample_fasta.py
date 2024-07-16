#!/usr/bin/env python3

import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq

def parse_arguments():
    parser = argparse.ArgumentParser(description="Simulate paired-end reads from a reference FASTA file.")
    parser.add_argument("reference", help="Path to the reference FASTA file")
    parser.add_argument("num_reads", type=int, help="Number of read pairs to simulate")
    parser.add_argument("read_length", type=int, help="Length of each read")
    parser.add_argument("insert_size", type=int, help="Average insert size between paired reads")
    parser.add_argument("output_prefix", help="Prefix for the output files")
    return parser.parse_args()

def read_fasta(reference_file):
    sequences = []
    for record in SeqIO.parse(reference_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def generate_paired_reads(sequences, num_reads, read_length, insert_size):
    paired_reads = []
    for _ in range(num_reads):
        seq = random.choice(sequences)
        start = random.randint(0, len(seq) - insert_size)
        forward_read = seq[start:start + read_length]
        reverse_read = Seq(seq[start + insert_size - read_length:start + insert_size]).reverse_complement()
        paired_reads.append((forward_read, reverse_read))
    return paired_reads

def write_paired_reads(paired_reads, output_prefix):
    forward_file = f"{output_prefix}_R1.fasta"
    reverse_file = f"{output_prefix}_R2.fasta"
    
    with open(forward_file, "w") as f1, open(reverse_file, "w") as f2:
        for i, (forward, reverse) in enumerate(paired_reads):
            f1.write(f">read_{i+1}/1\n{forward}\n")
            f2.write(f">read_{i+1}/2\n{reverse}\n")

def main():
    args = parse_arguments()
    sequences = read_fasta(args.reference)
    paired_reads = generate_paired_reads(sequences, args.num_reads, args.read_length, args.insert_size)
    write_paired_reads(paired_reads, args.output_prefix)

if __name__ == "__main__":
    main()
