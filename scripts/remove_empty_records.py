#!/usr/bin/env python
from Bio import SeqIO
from sys import argv
from pathlib import Path

path = Path(argv[1])
seq_records = []
for current_seq in SeqIO.parse(path, "fasta"):
    if len(current_seq.seq) != 0:
        seq_records.append(current_seq)

# Write the genome sequences to a file

new_path = path.parent / (path.name + ".tmp")

with open(new_path, "w") as output_handle:
    SeqIO.write(seq_records, output_handle, "fasta")
