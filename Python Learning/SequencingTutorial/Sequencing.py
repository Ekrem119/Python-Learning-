from Bio import SeqIO
import gzip
import matplotlib.pyplot as plt
import numpy as np
import os

file1 = "SRR800768_1_sub.fastq"
file2 = "SRR800768_2_sub.fastq"
## First Question
print("Counting reads")

def count_reads(filename):
    count = 0
    
    for record in SeqIO.parse(filename, "fastq"):
        count += 1
    return count


count1 = count_reads(file1)
print(f"{file1} has {count1} reads.")

count2 = count_reads(file2)
print(f"{file2} has {count2} reads.")

print("Verifying pairs")
records1 = SeqIO.parse(file1, "fastq")
records2 = SeqIO.parse(file2, "fastq")

for r1, r2 in zip(records1, records2):
    if r1.id.split('.')[0] != r2.id.split('.')[0]:
        print("Error: Mismatch!")
        break
else:
    print("match perfectly.")



## Second Question
filename = "SRR800768_1_sub.fastq"
first_record = next(SeqIO.parse(filename, "fastq"))

length = len(first_record.seq)
print(f"Read length is: {length} bp")

## Third Question



#records = SeqIO.parse('SRR800768_1_sub.fastq', 'fastq')
#for record in records:
 # seq = record.seq
  #phred_scores = record.letter_annotations['phred_quality']

  #print("Sequence:", seq) 