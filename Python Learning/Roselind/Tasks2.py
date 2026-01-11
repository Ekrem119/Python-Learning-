#task 1

s = "GCATTACGTCAGTACGTCAGTCAGTCAGTCAGTCAGTCAGCTAGACTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"

a_count =  s.count('A')
c_count =  s.count('C')
g_count =  s.count('G')
t_count =  s.count('T')
print(a_count,c_count, g_count, t_count)

#task1.2

s = "GCATTACGTCAGTACGTCAGTCAGTCAGTCAGTCAGTCAGCTAGACTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"

bases = ["A","C","G","T"]
counts = []

for base in bases:
    counts.append(str(s.count(base)))

print(counts)
