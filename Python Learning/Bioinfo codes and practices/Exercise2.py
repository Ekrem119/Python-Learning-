GENETIC_CODE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S",
    "TCA":"S","TCG":"S","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CTT":"L","CTC":"L",
    "CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R",
    "CGA":"R","CGG":"R","ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","AAT":"N","AAC":"N",
    "AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A",
    "GCA":"A","GCG":"A","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

raw_seq = "CAGCAGCTTGGCAACATTCTGTAAAAGCAGCTTTATACCTTTTAGCAAAGAAAAGGAGTTCGCCGGGCATAAAAGTAAAGGATGTCTTCTGGCAATTTCATATAAGTATTTTTTCAAAAATGTCTCTTCATTGTCATGAAAAGCAGTGCACATCATCAACCTCCTGGTCTCACCAATCGGGGGAGGTTTGGGTGTTCATCTTTGTTGTCAAGAAGGCATTTCATTTCTCTCAGGTTCTTGTTTTGCACAGCAGTCAGCCATTTCACCATAGGTTTCACGAAGAGTTGCAATGTGTCATAAATTTGTCTCCAAAAGGGTATGAAGTGATTTGTCACAATTTCAGCTGACTCATCAGCAACACATGTTTTTGCAAATTCACTTACTTCATTCACTAATTTTACATGATCTTCAAATGGACACTGCTGAAGACTGCAAGCAAGGCAATCAACACCAAGGCTTTGAAATTTTCTTCTCCCAAATCTTTAAACC GATGAGCAACCTCACTCTTGTGTCATGCTCGACGAAACACACCCCTGGAATAAGCCGAGCTAAAG"
dna = "".join(raw_seq.split()).upper()

all_protein_fragments = []
sequences_to_check = {
    "Forward": dna,
    "Reverse": dna.replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper()[::-1]
}

for name, seq in sequences_to_check.items():
    for frame in range(3):
        protein = ""
        for i in range(frame, len(seq), 3):
            codon = seq[i:i+3]
            if len(codon) == 3:
                protein += GENETIC_CODE.get(codon, "?")
        
        all_protein_fragments.extend(protein.split('*'))

longest_protein = max(all_protein_fragments, key=len)

print(f"Longest protein found: {len(longest_protein)} amino acids long.")
print(longest_protein)