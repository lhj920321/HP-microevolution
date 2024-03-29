
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser

def readRefRefGnm(RefRefGnmF):
	print("Loading Ref genome file :" + RefRefGnmF)
	RefGnm = ''
	for RefRefGnmFl in open(RefRefGnmF).readlines():
		if ">" in RefRefGnmFl:
			refID = RefRefGnmFl.split("\n")[0]
		else:
			RefGnm += RefRefGnmFl.split("\n")[0]

	return RefGnm

def CutGnm(RefGnm):
	out = []
	Qline = ''
	for x in range(0,cutLen):
		Qline += "F"
	for Posi in range(0,len(RefGnm),step):
		cutSeq = RefGnm[Posi:Posi + cutLen]	
		
		#out.append("@Posi_" + str(Posi+1) + "-" + str(Posi +  cutLen))
		#out.append(cutSeq)
		#out.append("+")
		#out.append(Qline)
		
		out.append(">Posi_" + str(Posi+1) + "-" + str(Posi +  cutLen))
		out.append(cutSeq)

	out = "\n".join(out)
	return out



def main():
	RefGnm = readRefRefGnm(RefRefGnmF)
	print("ref genome len : " + str(len(RefGnm)))
	out = CutGnm(RefGnm)
	#out_F = out_F + "/GnmeCutToReads.step" + str(step) + ".Len" + str(cutLen) + ".fasta"
	if (os.path.exists(out_F)):
		os.remove(out_F)
	out_SNP_F_O=open(out_F,'a')
	out_SNP_F_O.write(out + "\n")
	out_SNP_F_O.close()	






if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-r","--RefRefGnmF",
					  dest = "RefRefGnmF",
					  default = "",
					  metavar = "file",
					  help = "ref genome file .  [required]")

	parser.add_option("-o","--out_F",
					  dest = "out_F",
					  default = "",
					  metavar = "path",
					  help = "output stat file Path of snv Freq distribution.  [required]")
	parser.add_option("-s","--step",
					  dest = "step",
					  default = "",
					  metavar = "int",
					  help = "step to cut the ref Gnm.  [required]")

	parser.add_option("-l","--cutLen",
					  dest = "cutLen",
					  default = "",
					  metavar = "int",
					  help = "cutLen to cut the ref Gnm.  [required]")

	(options,args) = parser.parse_args()
	RefRefGnmF   = os.path.abspath(options.RefRefGnmF)
	out_F		  = os.path.abspath(options.out_F)
	step		  = int(options.step)
	cutLen		  = int(options.cutLen)


	#if ( not os.path.exists(out_F)):
		#os.mkdir(out_F)


	main()




