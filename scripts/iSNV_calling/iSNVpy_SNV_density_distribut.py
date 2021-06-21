
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 


def iSNV_SNP_tableRead(iSNV_SNP_table):
	iSNVSNPtablels = open(iSNV_SNP_table,'r').readlines()
	HeaderL = iSNVSNPtablels[0].split("\n")[0].split("\t")
	lie = 0
	lieD = {}
	for Header in HeaderL:
		lieD[Header] = lie
		lie += 1
	lsD = {}
	for iSNVSNPtablel in iSNVSNPtablels:
		if "#" not in iSNVSNPtablel and iSNVSNPtablel != "\n" :
			lL = iSNVSNPtablel.split("\n")[0].split("\t")
			Pos = lL[0]
			lsD[Pos] = lL
	return lieD,lsD


def samp_Freq_Dic(samples,lieD,lsD):
	samp_FreqDic = {}
	sampNA_PosiDic = {}
	#print(lieD)
	for sample in samples:
		Lie = lieD[sample]
		samp_FreqDic[sample] = {}
		sampNA_PosiDic[sample] = []
		for Posi in lsD:
			Freq = lsD[Posi][Lie]
			if Freq != "NA" :
				if Freq != "NO":				
					Freq = float(lsD[Posi][Lie])
					if Freq >= 0.02:
						samp_FreqDic[sample][Posi] = Freq
				else:
					Freq = 0
			else:
				sampNA_PosiDic[sample].append(Posi)
	#print(samp_FreqDic.keys())
	return samp_FreqDic,sampNA_PosiDic


def SNV_density(sample,samp_FreqDic,out_F):
	FreqDic = samp_FreqDic[sample]
	SNV_count_Dic = {}
	for Posi in FreqDic:
		Freq = FreqDic[Posi]
		if Freq >= minFreq and Freq < 1-minFreq:
			PosiFlag = int(int(Posi)/denst_denst_step)*denst_denst_step
			if PosiFlag not in SNV_count_Dic:
				SNV_count_Dic[PosiFlag] = 0
			SNV_count_Dic[PosiFlag] += 1
	print(SNV_count_Dic)
	densityLst,flagLst=[],[]
	outls=[]
	print(denst_denst_step)
	#for denstdenst_denst_step in SNV_count_Dic:		
	for denstdenst_denst_step in range(0,BacRefLen,denst_denst_step):		
		#print(denstdenst_denst_step)
		flagLst.append("P" + str(denstdenst_denst_step))
		if denstdenst_denst_step not in SNV_count_Dic :
			density = 0
			Count = 0
		else:
			density = round(float(SNV_count_Dic[denstdenst_denst_step])/denst_denst_step,2)
			Count = SNV_count_Dic[denstdenst_denst_step]
		densityLst.append(Count)
		outl = "\t".join([str(denstdenst_denst_step),str(Count)])
		outls.append(outl)

	print(len(flagLst))
	print("hhhh")

	outlines = "\n".join(outls) 
	out_F_O=open(out_F,'a')
	out_F_O.write(outlines)
	out_F_O.close()	
	
	out_plotF = out_F.split(".txt")[0] + ".pdf"
	plt.figure(dpi=300,figsize=(18,3))
	plt.bar(range(len(densityLst)), densityLst,color='red',tick_label=denstdenst_denst_step)  


	plt.title(sample.split(".sort")[0])
	plt.xlabel('Posi', fontsize = 13)
	plt.ylabel('count', fontsize = 13)
	plt.xticks(rotation=90)
	plt.savefig(out_plotF)
	plt.savefig(out_plotF.split(".pdf")[0] + ".png")







def main():

	lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
	#print(lieD)
	print("table : " + iSNV_SNP_table  +" iSNV table loaded!")
	iSNVtable_lst = list(lieD.keys())
	sampIDs = []
	for sampID in iSNVtable_lst:
		if "#Posi" not in sampID and "snv/SNP" not in sampID:
			sampIDs.append(sampID)

	print(str(len(sampIDs)) + "  samples are :"  + "\t".join(sampIDs))

	samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(sampIDs,lieD,lsD)

	for sample in sampIDs:
		out_F = out_P + "/" + sample + "iSNVpy.SNV_density_Stat" + str(minFreq) + "-" + str(1-minFreq) + ".txt"
		print("calculate for sample  : " + sample)
		print("output file :  " + out_F)
		snvSNPcount_Dic = SNV_density(sample,samp_FreqDic,out_F)





if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--iSNV_SNP_table",
					  dest = "iSNV_SNP_table",
					  default = "",
					  metavar = "file",
					  help = "iSNV table.  [required]")
	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output stat file Path of snv Freq distribution.  [required]")

	parser.add_option("-M","--minFreq",
					  dest = "minFreq",
					  default = "0",
					  metavar = "float",
					  help = "min Freq of snv to calculate (0-100%).  [required]")
	parser.add_option("-S","--denst_denst_step",
					  dest = "denst_denst_step",
					  default = "100",
					  metavar = "int",
					  help = "frame denst_denst_step to calculate.  [required]")
	parser.add_option("-L","--BacRefLen",
					  dest = "BacRefLen",
					  default = "1700000",
					  metavar = "int",
					  help = "Region of the bacteria genome to stat(default 1700000 for HP).  [required]")




	(options,args) = parser.parse_args()
	iSNV_SNP_table   = os.path.abspath(options.iSNV_SNP_table)
	out_P		  = os.path.abspath(options.out_P)
	minFreq		= float(options.minFreq)
	denst_denst_step	= int(options.denst_denst_step)
	BacRefLen   = int(options.BacRefLen)
	
	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




