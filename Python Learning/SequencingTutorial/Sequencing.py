from Bio import SeqIO
import gzip
import matplotlib.pyplot as plt
import numpy as np
import os

file1 = "SRR800768_1_sub.fastq"
file2 = "SRR800768_2_sub.fastq"

records = list(SeqIO.parse(file1, "fastq"))
read_length = len(records[0].seq)
num_reads = len(records)

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

base_counts = { 'A': np.zeros(read_length), 
                'T': np.zeros(read_length), 
                'C': np.zeros(read_length), 
                'G': np.zeros(read_length) }
for record in records:
    seq = str(record.seq)
    for i, base in enumerate(seq):
        if base in base_counts:
            base_counts[base][i] += 1
#Plot
plt.figure(figsize=(12, 6))
positions = np.arange(read_length)

colors = {'A': 'red', 'T': 'green', 'C': 'blue', 'G': 'black'}

for base, color in colors.items():
    percentages = (base_counts[base] / num_reads) * 100
    plt.plot(positions, percentages, color=color, label=base)

plt.xlabel("Position on Read")
plt.ylabel("Percentage %")
plt.title("Sequence Content across all bases")
plt.legend(loc='upper right')
plt.show()

## Fourth Question
total_bases = 0
gc_bases = 0

for record in records:
    seq = record.seq
    total_bases += len(seq)
    gc_bases += seq.count('G') + seq.count('C')

gc_content = (gc_bases / total_bases) * 100
print(f"Estimated GC Content: {gc_content:.2f}%")

## Fifth Question

qual_by_pos = [[] for _ in range(read_length)]

for record in records:
    phred_scores = record.letter_annotations["phred_quality"]
    for i, score in enumerate(phred_scores):
        if i < read_length: # Safety check
            qual_by_pos[i].append(score)

medians = []
q1s = []
q3s = []

for pos_scores in qual_by_pos:
    medians.append(np.median(pos_scores))
    q1s.append(np.percentile(pos_scores, 25))
    q3s.append(np.percentile(pos_scores, 75))

plt.figure(figsize=(12, 6))
plt.plot(positions, medians, color='purple', label='Median')
plt.plot(positions, q1s, color='blue', label='Q1 (25th)')
plt.plot(positions, q3s, color='red', label='Q3 (75th)')

plt.xlabel("Position")
plt.ylabel("Phred Quality Score")
plt.title("Quality Scores Distribution per Position")
plt.legend()
plt.grid(alpha=0.3)
plt.show()

""" def trim_read(record, quality_threshold=20):
    qualities = record.letter_annotations["phred_quality"]
    
    while qualities and qualities[-1] < quality_threshold:
        qualities.pop() # Remove last quality score
    
    return record[:len(qualities)]

f1_in = "SRR800768_1_sub.fastq"
f2_in = "SRR800768_2_sub.fastq"

out_p1 = open("paired_1.fastq", "w")
out_p2 = open("paired_2.fastq", "w")
out_u1 = open("unpaired_1.fastq", "w")
out_u2 = open("unpaired_2.fastq", "w")

min_len = 50
qual_thresh = 20

r1_iter = SeqIO.parse(f1_in, "fastq")
r2_iter = SeqIO.parse(f2_in, "fastq")

for r1, r2 in zip(r1_iter, r2_iter):
    
    r1_trimmed = trim_read(r1, qual_thresh)
    r2_trimmed = trim_read(r2, qual_thresh)
    
    r1_pass = len(r1_trimmed) >= min_len
    r2_pass = len(r2_trimmed) >= min_len
    
    if r1_pass and r2_pass:
        SeqIO.write(r1_trimmed, out_p1, "fastq")
        SeqIO.write(r2_trimmed, out_p2, "fastq")
        
    elif r1_pass and not r2_pass:
        SeqIO.write(r1_trimmed, out_u1, "fastq")
        
    elif not r1_pass and r2_pass:
        SeqIO.write(r2_trimmed, out_u2, "fastq")
        
out_p1.close()
out_p2.close()
out_u1.close()
out_u2.close()

print("Trimming and filtering complete.") """

#Plot for the cleaned data
print("Generating Quality Plot for Clean Data...")

clean_filename = "paired_1.fastq"
clean_records = list(SeqIO.parse(clean_filename, "fastq"))

max_len = max(len(r.seq) for r in clean_records)

clean_qual_by_pos = [[] for _ in range(max_len)]

for record in clean_records:
    scores = record.letter_annotations["phred_quality"]
    for i, score in enumerate(scores):
        clean_qual_by_pos[i].append(score)

c_medians = []
c_q1s = []
c_q3s = []
c_positions = []

for i, scores in enumerate(clean_qual_by_pos):
    if len(scores) > 0:
        c_positions.append(i)
        c_medians.append(np.median(scores))
        c_q1s.append(np.percentile(scores, 25))
        c_q3s.append(np.percentile(scores, 75))

# 5. Plotting
plt.figure(figsize=(12, 6))
plt.plot(c_positions, c_medians, color='purple', label='Median')
plt.plot(c_positions, c_q1s, color='blue', label='Q1 (25th)')
plt.plot(c_positions, c_q3s, color='red', label='Q3 (75th)')

plt.xlabel("Position")
plt.ylabel("Phred Quality Score")
plt.title("Quality Scores AFTER Cleaning (paired_1.fastq)")
plt.legend(loc='lower left')
plt.ylim(0, 45) 
plt.grid(alpha=0.3)
plt.show()