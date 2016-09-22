def cull(filenames):
	for filename in filenames:
		bases="ACGTN"
		shifts=[0]*16
		print shifts
		with open(filename, "rb") as f:
			
			with open(filename+"summary", "wb") as out:
				for line in f:
					line=line.split()
					if "nodes" in line:
						nodeschecked=float(line[3])	
					if "As" in line:
						acount=float(line[0])
						print acount
					elif "Cs" in line:
						ccount=float(line[0])
						print ccount
					elif "Gs" in line:
						gcount=float(line[0])
						print gcount
					elif "Ts" in line:
						tcount=float(line[0])
						print tcount

					if "->" in line:
						if line[0][1]!="N" and line[2]!="N":
							shifts[bases.index(line[0][1])*4+bases.index(line[2])]+=int(line[3])
							if int(line[3])>50:
								line[3]=str(int(line[3])/nodeschecked)
								out.write(" ".join(line)+"\n")
				out.write("\n")
				out.write("substitution frequency:")
				out.write(" "+str((sum(shifts)/nodeschecked)/85))#the way I got this number is that a read is 100bp. Barcode is 12. align is 5. 100-17=83	
				out.write("\n")
				for i in range(len(shifts)):
					if shifts[i]!=0:
						if i<3:
							print (bases[i/4]+" -> "+bases[i%4]+" "+str(shifts[i]/acount))
							out.write( bases[i/4]+" -> "+bases[i%4]+" "+str(shifts[i]/acount)+"\n")
						elif i<6:
							print (bases[i/4]+" -> "+bases[i%4]+" "+str(shifts[i]/ccount))
                                                        out.write( bases[i/4]+" -> "+bases[i%4]+" "+str(shifts[i]/ccount)+"\n")
						elif i<9:
							print (bases[i/4]+" -> "+bases[i%4]+" "+str(shifts[i]/gcount))
                                                        out.write( bases[i/4]+" -> "+bases[i%4]+" "+str(shifts[i]/gcount)+"\n")
						else:
							print (bases[i/4]+" -> "+bases[i%4]+" "+str(shifts[i]/tcount))
                                                        out.write( bases[i/4]+" -> "+bases[i%4]+" "+str(shifts[i]/tcount)+"\n")
cull(["all_0_thresh2.txt"])
					
