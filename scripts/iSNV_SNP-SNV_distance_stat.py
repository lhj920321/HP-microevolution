
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






def readVcf(VcfP,persNAPosiD):
	files = os.listdir(VcfP)
	GeneSNPCountDic = {}
	SNP_distanceL,SNV_distanceL = [],[]

	sample_Stat = {}
	sampTypeCount = {}
	for File in files:
		SNP_Cites,SNV_Cites = [],[]

		if "P" in File:	
			VcfF = VcfP + "/" + File
			PersID = File.split("-")[0]
			PersNAPosi = persNAPosiD[PersID]
			SampID = File.split(".minFreq0.02.snpEffAnno.vcf")[0].split("_2_")[0]		
			#refID = File.split(".minFreq0.02.snpEffAnno.vcf")[0].split("_2_")[1]
			print("calculate for:  " + SampID)
			VcfFls = open(VcfF).readlines()
			sample_Stat[SampID] = []
			sampTypeCount[SampID] = {}
			for VcfFl in VcfFls:
				if "#" not in VcfFl:
					Tags = VcfFl.split("\n")[0].split("\t")
					Cite = int(Tags[1])
					Freq = float(Tags[9].split(":")[5].split("%")[0])
					#print(Freq)
					if Cite not in PersNAPosi:
						if Freq >= minFreqofSNP * 100:
							if "SNP" not in sampTypeCount[SampID]:
								sampTypeCount[SampID]["SNP"] = 0
							sampTypeCount[SampID]["SNP"] += 1
							if Cite not in SNP_Cites:
								SNP_Cites.append(Cite)
							if len(SNP_Cites) != 1:
								SNP_Dicstans = SNP_Cites[-1] - SNP_Cites[-2]
								if SNP_Dicstans < 1000000:
									SNP_distanceL.append("\t".join([SampID,"SNP",str(SNP_Dicstans)]))
									sample_Stat[SampID].append("\t".join([SampID,"SNP",str(SNP_Dicstans)]))
						elif Freq >= (1-minFreqofSNP) * 100 and Freq < minFreqofSNP * 100:
							if "SNV" not in sampTypeCount[SampID]:
								sampTypeCount[SampID]["SNV"] = 0
							sampTypeCount[SampID]["SNV"] += 1
							if Cite not in SNV_Cites:
								SNV_Cites.append(Cite)

							if len(SNV_Cites) != 1:
								#print(SNV_Cites[-1],SNV_Cites[-2])
								#print()
								SNV_Dicstans = SNV_Cites[-1] - SNV_Cites[-2]
								if SNV_Dicstans < 1000000:
									SNV_distanceL.append("\t".join([SampID,"SNV",str(SNV_Dicstans)]))
									sample_Stat[SampID].append("\t".join([SampID,"SNV",str(SNV_Dicstans)]))
	outs = ["\t".join(["sampID","type","Distance"])]
	outs.append("\n".join(SNP_distanceL))
	outs.append("\n".join(SNV_distanceL))
	OUT = "\n".join(outs)

	print(OUT)


	outF = out_P + "/SNPminFreq"+ str(minFreqofSNP) + "_distance.txt"
	if (os.path.exists(outF)):
		os.remove(outF)
	out_P_O=open(outF,'a')
	out_P_O.write(OUT + "\n")
	out_P_O.close()	


	print(sampTypeCount)
	for SampID in sample_Stat:
		if len(sampTypeCount[SampID]) ==2 and sampTypeCount[SampID]["SNP"] >= 20 and sampTypeCount[SampID]["SNV"] >= 20:
			outs = ["\t".join(["sampID","type","Distance"])]
			outs.append("\n".join(sample_Stat[SampID]))
			OUT = "\n".join(outs)

			outF = out_P + "/SNPminFreq"+ str(minFreqofSNP) + "." + SampID + ".P_distance.txt"
			if (os.path.exists(outF)):
				os.remove(outF)
			out_P_O=open(outF,'a')
			out_P_O.write(OUT + "\n")
			out_P_O.close()	





def main():
	persNAPosiD = {}
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
			persNAPosiD[personID] = personNAPosi
	print(vcfP)
	readVcf(vcfP,persNAPosiD)
	

	#out_P_O.write("\n".join(OUts) + "\n")
	#out_P_O.close()	





if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--iSNV_SNP_tableP",
					  dest = "iSNV_SNP_tableP",
					  default = "",
					  metavar = "file",
					  help = "iSNV table.  [required]")
	parser.add_option("-v","--vcfP",
					  dest = "vcfP",
					  default = "",
					  metavar = "path",
					  help = "Path of snpeff vcf files .  [required]")


	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "file",
					  help = "output stat file Path of snv Freq distribution.  [required]")

	parser.add_option("-m","--minFreqofSNP",
					  dest = "minFreqofSNP",
					  default = "0.95",
					  metavar = "float",
					  help = "min Freq of SNP to calculate (0-1).  [required]")


	(options,args) = parser.parse_args()
	iSNV_SNP_tableP   = os.path.abspath(options.iSNV_SNP_tableP)
	vcfP    = os.path.abspath(options.vcfP)
	out_P		  = os.path.abspath(options.out_P)
	minFreqofSNP		= float(options.minFreqofSNP)


	print(vcfP)
	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()





