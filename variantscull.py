def cull(filenames):
	for filename in filenames:
		bases="ACGTN"
		shifts=[0]*16
		print shifts
		with open(filename, "rb") as f:
			with open(filename+"summary", "wb") as out:
				for line in f:
					if "->" in line:
						line=line.split()
						if line[0][1]!="N" and line[2]!="N":
							shifts[bases.index(line[0][1])*4+bases.index(line[2])]+=int(line[3])
							if int(line[3])>50:
								line[3]=str(float(line[3])/377183)	
								out.write(" ".join(line)+"\n")
				out.write("\n")
				for i in range(len(shifts)):
					if shifts[i]!=0:
						out.write( bases[i/4]+" -> "+bases[i%4]+" "+str(shifts[i]/377183.0)+"\n")


cull(["all_0_thresh7.txt"])
					
