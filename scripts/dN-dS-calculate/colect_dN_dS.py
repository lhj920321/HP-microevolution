
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 



def codeml_run(codemlP,Hl,FiltSamps):
	codemlFs = []
	for file in os.listdir(codemlP):
		if ".codeml" in file:
			codemlFs.append(codemlP + "/" + file)
	#print(codemlFs)

	OUTs,OUTPositive = [],[]
	OUTs.append(Hl)
	OUTPositive.append(Hl)

	for codemlF in codemlFs:
		codemlFls = open(codemlF).readlines()
		if codemlFls != []:
			lineC,PositiveC = 0,0
			lineD = {}
			fileLs = []
			for codemlFl in codemlFls:
				lineD[lineC] = codemlFl
				if "dN/dS=" in codemlFl:
					Tags = codemlFl.split("\n")[0].split("=")
					print(Tags)
					S = str(float(Tags[2].split("N")[0]))
					N = str(float(Tags[3].split("dN/dS")[0]))
					dNdS = str(float(Tags[4].split("dN")[0]))
					dN = str(float(Tags[5].split("dS")[0]))
					dS = str(float(Tags[6]))

					if float(dNdS) >= 10:
						dNdS="0"
					Sampsl = lineD[lineC-4]
					Samp1 = Sampsl.split("...")[0].split("(")[1].split(")")[0]
					Samp2 = Sampsl.split("...")[1].split("(")[1].split(")")[0]
					GeneGrp =  Samp1.split("__")[0]
					GeneID1 =  Samp1.split("__")[2]
					GeneID2 =  Samp2.split("__")[2]
					SampID0 = Samp1.split("__")[1].split("-")[0] 
					SampGrp = Samp1.split("__")[1].split("-")[1].split("_")[0]  

					SampID1 = Samp1.split("__")[1].split("-")[2]
					SampID2 = Samp2.split("__")[1].split("-")[2]
#dN > 0, dS < 2 and dN/dS < 10
					if Samp1.split("__")[1] not in FiltSamps and (Samp2.split("__")[1] not in FiltSamps) :
						SampIDs = [SampID1,SampID2]
						SampIDs.sort()
						print(SampIDs)
						print("_".join(SampIDs))
						#if float(dN) >0  and float(dS) <2 and float(dNdS) < 50 and  float(dNdS) != 0:
						personFM = str("%02d" % int(SampID0.split("P")[1]))
						if float(dNdS) < 50  and float(dS) < 2 :  #and  float(dNdS) != 0
							PositiveC += 1
							outl = "\t".join(["P" + personFM,SampGrp,SampID0 + "-" + "_".join(SampIDs),GeneGrp,GeneID1 + "__" + GeneID2,S,N,dNdS,dN,dS])
							OUTs.append(outl)
							fileLs.append(outl)



				lineC += 1

			if PositiveC >= 1:
				OUTPositive = OUTPositive + fileLs
			
	OUT = "\n".join(OUTs)
	OUTposit = "\n".join(OUTPositive)
	return OUT,OUTposit




					





def main():
	FiltSamps=["P1-E-x","P2-E-x","P3-E-t","P4-E-x","P6-E-d","P7-E-6","P8-E-d","P8-E-x","P9-E-d","P10-E-x","P10-E-x","P23-C_sL-x","P24-C_sC-x"]
	FiltSamps=[]
	Hl = "\t".join(["person","Samp_Group","Samp1__Samp2","GeneGrp","Samp1Gene__Samp1Gene","S","N","dNdS","dN","dS"])

	OUT,OUTposit = codeml_run(codemlP,Hl,FiltSamps)

	out_F_O=open(out_F,'a')
	out_F_O.write(OUT + "\n")
	out_F_O.close()		



'''
	F = out_F.split(".txt")[0] + ".over1.txt"
	if ( os.path.exists(F)):
		os.remove(F)
	out_F_O=open(F,'a')
	print(OUTposit)
	out_F_O.write(OUTposit + "\n")
	out_F_O.close()		
'''




if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)


	parser.add_option("-c","--codemlP",
					  dest = "codemlP",
					  default = "",
					  metavar = "path",
					  help = "Path of codeml files .  [required]")


	parser.add_option("-o","--out_F",
					  dest = "out_F",
					  default = "",
					  metavar = "path",
					  help = "output stat file Path of snv Freq distribution.  [required]")



	(options,args) = parser.parse_args()
	codemlP   = os.path.abspath(options.codemlP)
	out_F		  = os.path.abspath(options.out_F)



	if ( os.path.exists(out_F)):
		os.remove(out_F)



	main()





