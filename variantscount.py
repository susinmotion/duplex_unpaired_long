variants=[0]*256
bases="ACGT"
with open("output2.txt") as f:
	for line in f:
		try:
			variants[int(line)]+=1
		except:
			continue

for i in range(len(variants)):
	base1=bases[i/64]
	base2=bases[i%64/16]
	base3=bases[i%16/4]
	shift=bases[i%4]
	trio=base1+base2+base3
	print trio, "->", shift, variants[i]
