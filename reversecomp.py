def ReverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def makeNewFile(file):
    with open(file, "rb") as f:
	with open(str(file)+"rev", "wb") as ff:
	   line=f.readline()
	   while(line):
		ff.write(line)
		line=f.readline()
		print line
		ff.write(ReverseComplement(line.strip("\n").strip())+"\n")
		ff.write(f.readline())
		ff.write(f.readline())
		line=f.readline()

makeNewFile("veryshort5.txt")
