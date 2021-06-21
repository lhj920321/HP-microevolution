
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
					#if Freq >= 0.02:
					samp_FreqDic[sample][Posi] = Freq
				else:
					Freq = 0
					samp_FreqDic[sample][Posi] = Freq
			else:
				sampNA_PosiDic[sample].append(Posi)
	return samp_FreqDic,sampNA_PosiDic



def Person_NA_Posi(samples,sampNA_PosiDic):
	personNAPosi = {}
	for sample in samples:
		sampleNAPosi = sampNA_PosiDic[sample]
		for posi in sampleNAPosi:
			if posi not in personNAPosi:
				personNAPosi[posi] = ''
	return personNAPosi



def Person_SNP_with_SNV_Posi(samples,samp_FreqDic,lsD):
	personSNP_with_SNVCite = {}
	for posi in lsD.keys():
		count = 0
		for sample in samples:
			sampleFreqDic = samp_FreqDic[sample]
			if posi in sampleFreqDic:
				FREQ = sampleFreqDic[posi]
				if FREQ >= 0.02 and FREQ < SNPminFreq:
					count += 1
		if count >= 1:
			personSNP_with_SNVCite[posi] = ''
			
	return personSNP_with_SNVCite



def Compare_diffSNPLOci(samp1,samp2,samp_FreqDic,personNAPosi,lsD,subgroup,group):
	count = 0
	refSamp = samp1.split("_2_")[1]
	for Posi in lsD:
		if Posi not in personNAPosi:
			Samp1Frq = samp_FreqDic[samp1][Posi]
			Samp2Frq = samp_FreqDic[samp2][Posi]
			if (Samp1Frq >= SNPminFreq and Samp2Frq < 0.02) or (Samp2Frq >= SNPminFreq and Samp1Frq < 0.02):
				count += 1

	return "\t".join([refSamp,samp1.split("_2_")[0] + "--" + samp2.split("_2_")[0],str(count),subgroup,group])




def Compare_diffSNPLOci_noSNV(samp1,samp2,samp_FreqDic,personNAPosi,personSNP_with_SNVCite,lsD,subgroup,group):
	count = 0
	refSamp = samp1.split("_2_")[1]
	for Posi in lsD:
		if Posi not in personNAPosi and Posi not in personSNP_with_SNVCite:
			Samp1Frq = samp_FreqDic[samp1][Posi]
			Samp2Frq = samp_FreqDic[samp2][Posi]
			if (Samp1Frq >= SNPminFreq and Samp2Frq < 0.02) or (Samp2Frq >= SNPminFreq and Samp1Frq <0.02):
				count += 1

			
	return "\t".join([refSamp,samp1.split("_2_")[0] + "--" + samp2.split("_2_")[0],str(count),subgroup,group])





def main():
	lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
	print(lieD)
	iSNVtable_lst = list(lieD.keys())
	sampIDs = []
	for sampID in iSNVtable_lst:
		if "#Posi" not in sampID and "snv/SNP" not in sampID and sampID not in ["P24-C_sC-j_2_P24-C_sC-t","P8-E-t_2_P8-E-x"]:
			sampIDs.append(sampID)
	print(sampIDs)
	samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(sampIDs,lieD,lsD)
	personNAPosi = Person_NA_Posi(sampIDs,sampNA_PosiDic)

	personSNP_with_SNVCite = Person_SNP_with_SNV_Posi(sampIDs,samp_FreqDic,lsD)

	out_F_O=open(out_F,'a')
	for Idx in range(0,len(sampIDs)):
		Samp2Idx = Idx + 1
		Samp1=sampIDs[Idx]	
		Samp2Lst = sampIDs[Samp2Idx:len(sampIDs)]
		subgroup = Samp1.split("-")[1]
		if "E" in subgroup:
			group = subgroup 
		else:
			group = subgroup.split("_")[0] 

		for Samp2 in Samp2Lst:
			#stat = Compare_diffSNPLOci(Samp1,Samp2,samp_FreqDic,personNAPosi,lsD,subgroup,group)
			##without SNV(0.02-0.95)
			stat = Compare_diffSNPLOci_noSNV(Samp1,Samp2,samp_FreqDic,personNAPosi,personSNP_with_SNVCite,lsD,subgroup,group)


			out_F_O.write(stat + "\n")
	out_F_O.close()	


	
'''
	personCommonSNPCite = Person_commonSNP_Posi(sampIDs,samp_FreqDic,lsD)
	snvSNPcount_Dic = {}

	for sample in sampIDs:
		#print(sample)
		out_F = out_P + "/" + sample + "iSNVpy.Freq_distribu.txt"
		snvSNPcount_Dic = tongjiFREQ(sample,samp_FreqDic,out_F,personNAPosi,personCommonSNPCite,snvSNPcount_Dic)

	person = iSNV_SNP_table.split("/")[-2]
	out_snvSNPCountF = out_P + "/" + person + ".snvSNPcount.txt"
	out_snvSNPCountF_O=open(out_snvSNPCountF,'a')
	out_snvSNPCountF_O.write("\t".join(["sample","snvCount(" + str(minFreq) + "-" + str(1-minFreq) + ")" ,"SNPCount(>=" + str(1-minFreq) + ")","personCommonCitesNum"]) + "\n")
	for sampleID in snvSNPcount_Dic:	
		out_snvSNPCountF_O.write(snvSNPcount_Dic[sampleID] + "\n")
	out_snvSNPCountF_O.close()	

'''




if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--iSNV_SNP_table",
					  dest = "iSNV_SNP_table",
					  default = "",
					  metavar = "file",
					  help = "iSNV table.  [required]")
	parser.add_option("-o","--out_F",
					  dest = "out_F",
					  default = "",
					  metavar = "file",
					  help = "output stat file Path of snv Freq distribution.  [required]")

	parser.add_option("-m","--SNPminFreq",
					  dest = "SNPminFreq",
					  default = "0.98",
					  metavar = "float",
					  help = "min Freq of SNP to calculate (0-1).  [required]")




	(options,args) = parser.parse_args()
	iSNV_SNP_table   = os.path.abspath(options.iSNV_SNP_table)
	out_F		  = os.path.abspath(options.out_F)
	SNPminFreq		= float(options.SNPminFreq)



	if ( os.path.exists(out_F)):
		#os.mkdir(out_F)
		os.remove(out_F)


	main()





