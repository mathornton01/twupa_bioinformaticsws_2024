import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq

def introduce_errors(seq, ins_prob, del_prob, mut_prob):
    bases = ['A', 'C', 'G', 'T']
    new_seq = []
    i = 0

    while i < len(seq):
        r = random.random()
        if r < ins_prob:
            # Insertion
            new_seq.append(random.choice(bases))
        elif r < ins_prob + del_prob:
            # Deletion
            i += 1
            continue
        elif r < ins_prob + del_prob + mut_prob:
            # Mutation
            new_seq.append(random.choice([b for b in bases if b != seq[i]]))
        else:
            new_seq.append(seq[i])
        i += 1

    return ''.join(new_seq)

def simulate_paired_end_reads(reference, num_reads, read_length, insert_length, ins_prob, del_prob, mut_prob):
    records = list(SeqIO.parse(reference, "fasta"))
    paired_reads = []

    for _ in range(num_reads):
        record = random.choice(records)
        seq_len = len(record.seq)

        if seq_len < insert_length:
            continue

        start_pos = random.randint(0, seq_len - insert_length)
        insert_seq = record.seq[start_pos:start_pos + insert_length]

        insert_seq = introduce_errors(str(insert_seq), ins_prob, del_prob, mut_prob)

        read1 = Seq(insert_seq[:read_length])
        read2 = Seq(insert_seq[-read_length:]).reverse_complement()

        paired_reads.append((read1, read2))

    return paired_reads

def write_fastq(paired_reads, output_prefix):
    with open(f"{output_prefix}_1.fastq", "w") as f1, open(f"{output_prefix}_2.fastq", "w") as f2:
        for i, (read1, read2) in enumerate(paired_reads):
            f1.write(f"@read{i}/1\n{read1}\n+\n{'~' * len(read1)}\n")
            f2.write(f"@read{i}/2\n{read2}\n+\n{'~' * len(read2)}\n")

def main():
    parser = argparse.ArgumentParser(description="Simulate paired-end reads from a reference FASTA file with errors.")
    parser.add_argument("reference", type=str, help="Path to the reference FASTA file")
    parser.add_argument("num_reads", type=int, help="Number of reads to simulate")
    parser.add_argument("read_length", type=int, help="Length of each read")
    parser.add_argument("insert_length", type=int, help="Length of the insert")
    parser.add_argument("output_prefix", type=str, help="Prefix for the output FASTQ files")
    parser.add_argument("--ins_prob", type=float, default=0.001, help="Probability of insertion errors")
    parser.add_argument("--del_prob", type=float, default=0.001, help="Probability of deletion errors")
    parser.add_argument("--mut_prob", type=float, default=0.001, help="Probability of mutation errors")

    args = parser.parse_args()

    paired_reads = simulate_paired_end_reads(args.reference, args.num_reads, args.read_length, args.insert_length, args.ins_prob, args.del_prob, args.mut_prob)
    write_fastq(paired_reads, args.output_prefix)

if __name__ == "__main__":
    main()

