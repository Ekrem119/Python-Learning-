""" #task1 
import this

#task2
a = 799
b = 800
c = (a*a)+(b*b)
#alternative task 2.1 
a = 799
b = 800
c = (a**2 + b**2)
#alternative task 2.2
with open('rosalind_ini2.txt' , 'r') as f:
    line = f.read()
a, b = map(int, line.split())
print(a**2 + b**2)    

#task 3 
s = "gKt5gTflBEyeRF9JBBcjNRIAlausrBp4f2pzVBaZp8L9P6lrXJolIPfd0XaiXEk8sWPvjcd8XHG00NtaealBRRmykissKnnKPzXdEQGS8jazRw4pLamVG4bpAg4wSoByHZ1ko6jrRUM4VO6gorSjcax6nqaSoJtQyXkbRSOV6VgdCTE1G7iGcPMJThjOYTou"
slice = s[23:27+1]  
slice2 = s[86:91+1]
print(slice + " " + slice2) """

""" #task 4
a = 4633
b = 8808

sum_of_odds = 0

for i in range(a, b+1):
    if i % 2 == 1:
        sum_of_odds += i

print(sum_of_odds) """

""" #task 5

with open('rosalind_ini5.txt', 'r') as f:

    for line_number, line in enumerate(f,1):
        if line_number %2 == 0:
            print(line.strip())
        
 """
#task 6
s = "When I find myself in times of trouble Mother Mary comes to me Speaking words of wisdom let it be And in my hour of darkness she is standing right in front of me Speaking words of wisdom let it be Let it be let it be let it be let it be Whisper words of wisdom let it be And when the broken hearted people living in the world agree There will be an answer let it be For though they may be parted there is still a chance that they will see There will be an answer let it be Let it be let it be let it be let it be There will be an answer let it be Let it be let it be let it be let it be Whisper words of wisdom let it be Let it be let it be let it be let it be Whisper words of wisdom let it be And when the night is cloudy there is still a light that shines on me Shine until tomorrow let it be I wake up to the sound of music Mother Mary comes to me Speaking words of wisdom let it be Let it be let it be let it be yeah let it be There will be an answer let it be Let it be let it be let it be yeah let it be Whisper words of wisdom let it be"

split_s = s.split()

d = {}

for word in split_s:
    d[word] = d.get(word, 0) + 1

for word, count in d.items():
    print(f"{word} {count}")