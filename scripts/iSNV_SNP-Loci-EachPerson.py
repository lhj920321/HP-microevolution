
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



def Person_SNPLoci(personID,group,subgroup,samples,samp_FreqDic,lsD,personNAPosi):
	personSNPLoci = {}
	for posi in lsD.keys():
		count = 0
		for sample in samples:
			sampleFreqDic = samp_FreqDic[sample]
			if posi not in personNAPosi and posi in sampleFreqDic:
				FREQ = sampleFreqDic[posi]
				if FREQ >= SNPminFreq and posi not in personSNPLoci:
					personSNPLoci[posi] = ''
				
			
	return "\t".join([personID,str(len(personSNPLoci)),group,subgroup])




def main():
	OUts =[]
	OUts.append("\t".join(["person","SNPLoci","group","subgroup"]))
	out_F_O=open(out_F,'a')

	for Person in range(1,26):
		if Person not in [13,14]:
			iSNV_SNP_table = iSNV_SNP_tableP + "/P" + str(Person)  + "/all.iSNV_with_SNP.pyResults.txt" 
			lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
			
			iSNVtable_lst = list(lieD.keys())
			sampIDs = []
			for sampID in iSNVtable_lst:
				if "#Posi" not in sampID and "snv/SNP" not in sampID and sampID not in ["P24-C_sC-j_2_P24-C_sC-t","P8-E-t_2_P8-E-x"]:
					sampIDs.append(sampID)
			print(Person)
			print(sampIDs)
			samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(sampIDs,lieD,lsD)
			personNAPosi = Person_NA_Posi(sampIDs,sampNA_PosiDic)

			subgroup = sampIDs[0].split("-")[1]
			personID = sampIDs[0].split("-")[0]
			if "E" in subgroup:
				group = subgroup 
			else:
				group = subgroup.split("_")[0] 

			personID = "P" + str(Person)
			stat = Person_SNPLoci(personID,group,subgroup,sampIDs,samp_FreqDic,lsD,personNAPosi)
			OUts.append(stat)

	out_F_O.write("\n".join(OUts) + "\n")
	out_F_O.close()	





if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--iSNV_SNP_tableP",
					  dest = "iSNV_SNP_tableP",
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
					  default = "0.9",
					  metavar = "float",
					  help = "min Freq of SNP to calculate (0-1).  [required]")


	(options,args) = parser.parse_args()
	iSNV_SNP_tableP   = os.path.abspath(options.iSNV_SNP_tableP)
	out_F		  = os.path.abspath(options.out_F)
	SNPminFreq		= float(options.SNPminFreq)



	if ( os.path.exists(out_F)):
		os.remove(out_F)


	main()





