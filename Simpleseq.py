#generate n mers

import random

n = 50
''.join(random.choice('AGTC') for i in range(n))

#count number of specific pattern in a text

seq = 'CGACCAGACTGCAATGTTCGC'
seq.count('ATG') #count 'ATG' string in a sequence

#GC content

GC = (seq.count('G')+seq.count('C'))*100/len(seq) 

#find all position of pattern in string

def findpos(pattern, sequence):
	pos = []
	i = sequence.find(pattern)
	while i != -1:
		pos.append(i)
		i = sequence.find(pattern, i+1)
	for j in range(len(pos)):
		pos[j] += 1
	print ('The pattern,', pattern, ',can be aligned at', sequence.count(pattern), 'positions:', pos)

#e.g. findall('ATGTGG', seq)