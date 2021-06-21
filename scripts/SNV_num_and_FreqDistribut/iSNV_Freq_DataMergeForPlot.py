
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 


def readF(F):
	Fls = open(F,'r').readlines()
	ls = []
	for Fl in Fls:
		if Fl != "\n" :
			Tgs = Fl.split("\n")[0].split("\t")
			#print("\t".join([Tgs.split("_2_")[0],Tgs[1],Tgs[2]]))
			ls.append("\t".join([Tgs[0].split("_2_")[0],Tgs[1],Tgs[2],"yes"]))
			#ls.append(Fl)
	Ols = "\n".join(ls)
	return Ols

def outputNullLines(buweiFileName):
	outls = []
	for Idx in range(0,101,2):
		outl = "\t".join([buweiFileName.split("_2_")[0],str(Idx),"0","no"])
		outls.append(outl)
	NULLls = "\n".join(outls)
	return NULLls


def main():
	buweiLst = ["x","d","j","t"]
	for per in range(1,26):
		if per not in [13,14]:
			Person = "P" + str(per)
			persPath = NGS_iSNV_stat_P + "/" + Person + "/SNVFreq_distribut_Stat_filtRepeat/"
			Files = os.listdir(persPath)
			SampFiles = []
			for File in Files:
				if ".Freq_distribu.txt" in File and File.split(".iSNVpy")[1] == ".Freq_distribu.txt" :
					SampFiles.append(File)	
			print(SampFiles)
			Tags = SampFiles[0].split("_2_")[0].split("-")
			out_F = out_P + "/" + Person + ".iSNV_0.95.Freq_distrib.txt"
			if (os.path.exists(out_F)):
				os.remove(out_F)
			print(SampFiles)
			out_F_O = open(out_F,'a')
			for buwei in buweiLst:
				print(buwei)
				buweiFileName = "-".join(Tags[0:2] + [buwei]) +  "_2_" + SampFiles[0].split("_2_")[1]
				print(buweiFileName)
				if str(buweiFileName) in SampFiles:
					F = persPath  + buweiFileName
					Ls = readF(F)
				else:
					Ls = outputNullLines(buweiFileName)
					print("nono")
				#print(Ls)
			
				
				out_F_O.write(Ls + "\n")
			out_F_O.close()	




if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--NGS_iSNV_stat_P",
					  dest = "NGS_iSNV_stat_P",
					  default = "",
					  metavar = "path",
					  help = "NGS_iSNV_stat_P.  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output path .  [required]")


	(options,args) = parser.parse_args()
	NGS_iSNV_stat_P   = os.path.abspath(options.NGS_iSNV_stat_P)
	out_P		  = os.path.abspath(options.out_P)



	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()





