def get_reverse_complement(dna_sequence):
    s = dna_sequence.upper()
    
    if any(c not in 'ATCG' for c in s):
        raise ValueError("Invalid DNA sequence. Must only contain A, T, C, or G.")
    
    comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    complement = "".join(comp_map[c] for c in s)
    
    return complement[::-1]

if __name__ == "__main__":
    try:
        dna_input = input("Enter DNA sequence: ")
        result = get_reverse_complement(dna_input)
        print(f"Reverse Complement: {result}")
    except ValueError as e:
        print(f"Error: {e}")