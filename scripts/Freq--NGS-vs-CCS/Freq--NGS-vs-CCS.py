#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys,os
from optparse import OptionParser



def iSNV_SNP_tableRead(iSNV_SNP_table):
	print(iSNV_SNP_table)
	if os.path.exists(iSNV_SNP_table):
		print(iSNV_SNP_table)
	iSNVSNPtablels = open(iSNV_SNP_table).readlines()
	HeaderL = iSNVSNPtablels[0].split("\n")[0].split("\t")
	lie = 0
	lieD = {}
	for Header in HeaderL:
		if "P8-E-t" not in Header and ("P24-C_sC-j" not in Header) and ("P11-E-j" not in Header):
			lieD[Header] = lie
		lie += 1
	lsD = {}
	for iSNVSNPtablel in iSNVSNPtablels:
		if "#" not in iSNVSNPtablel and iSNVSNPtablel != "\n" :
			lL = iSNVSNPtablel.split("\n")[0].split("\t")
			Pos = lL[0]
			lsD[Pos] = lL
	Samps = []
	for S in lieD:
		if S != "#Posi" and S != "snv/SNP":
			Samps.append(S)
	print(lieD)
	P_S_FreqD = {}
	for SP in Samps:
		P_S_FreqD[SP] = {}
		Lie = lieD[SP]
		for Posi in lsD:
			Freq = lsD[Posi][Lie]
			if Freq != "NA" :
				if Freq != "NO":				
					Freq = lsD[Posi][Lie]
				else:
					Freq = "0"
					
				if float(Freq) != 0:
					P_S_FreqD[SP][Posi] = Freq
			else:
				P_S_FreqD[SP][Posi] = Freq

	return P_S_FreqD


def main():
	CCS_Samps = ["P1-E-d_2_P1-E-j",\
	"P1-E-j_2_P1-E-j","P1-E-t_2_P1-E-j","P1-E-x_2_P1-E-j",\
	"P2-E-t_2_P2-E-t",\
	"P4-E-d_2_P4-E-j","P4-E-j_2_P4-E-j","P4-E-x_2_P4-E-j",\
	"P5-E-d_2_P5-E-j","P6-E-d_2_P6-E-t","P6-E-t_2_P6-E-t",\
	"P7-E-d_2_P7-E-j","P7-E-j_2_P7-E-j""P18-C_bR-j_2_P18-C_bR-x",\
	"P18-C_bR-x_2_P18-C_bR-x","P19-C_bR-d_2_P19-C_bR-j",\
	"P19-C_bR-j_2_P19-C_bR-j","P19-C_bR-t_2_P19-C_bR-j",\
		"P20-C_bR-d_2_P20-C_bR-t","P20-C_bR-j_2_P20-C_bR-t",\
		"P20-C_bR-t_2_P20-C_bR-t","P20-C_bR-x_2_P20-C_bR-t",\
		"P21-C_sL-d_2_P21-C_sL-j","P21-C_sL-j_2_P21-C_sL-j",\
		"P21-C_sL-t_2_P21-C_sL-j","P21-C_sL-x_2_P21-C_sL-j",\
		"P22-C_sL-d_2_P22-C_sL-d","P22-C_sL-t_2_P22-C_sL-d",\
		"P23-C_sL-j_2_P23-C_sL-j","P23-C_sL-x_2_P23-C_sL-j",\
		"P24-C_sC-j_2_P24-C_sC-t","P24-C_sC-t_2_P24-C_sC-t",\
		"P24-C_sC-x_2_P24-C_sC-t","P25-C_sC-j_2_P25-C_sC-t",\
		"P25-C_sC-t_2_P25-C_sC-t",]

	CCSpersons = ["P1","P2","P4","P5","P6","P7","P8","P18","P19","P20","P21","P22","P23","P24","P25"]
	#CCSpersons = ["P1"]

	print(len(CCS_Samps))

	out_F = out_P + "/All-CCS-vs-NGS.Freq.txt"
	if (os.path.exists(out_F)):
		os.remove(out_F)
	out_FO = open(out_F,'a')
	for CCSperson in CCSpersons:
		print(CCSperson)
		CCSTable = CCS_table_P + "/" + CCSperson + "/all.iSNV_with_SNP.pyResults.txt"
		print(CCSTable)
		print(CCS_table_P)
		CCS_P_S_FreqD = iSNV_SNP_tableRead(CCSTable)
		NGSTable = NGS_table_P + "/" + CCSperson + "/all.iSNV_with_SNP.pyResults.txt"
		print(NGSTable)
		NGS_P_S_FreqD = iSNV_SNP_tableRead(NGSTable)
		
		for CCS_Samp in CCS_Samps:
			if CCS_Samp.split("-")[0] == CCSperson:
				allPosi = {}
				Os = []
				for CCS_Po in CCS_P_S_FreqD[CCS_Samp]:
					if CCS_Po not in allPosi:
						allPosi[CCS_Po] = ''
				for ngs_Po in NGS_P_S_FreqD[CCS_Samp] :
					if ngs_Po not in allPosi:
						allPosi[ngs_Po] = ''

				for Po in allPosi:
					countNA,count0CCS,count0NGS = 0,0,0
					ccs_Freq,ngs_Freq = "0","0"
					if Po in CCS_P_S_FreqD[CCS_Samp]:
						ccs_Freq = CCS_P_S_FreqD[CCS_Samp][Po]
						if ccs_Freq == "NA":
							countNA += 1
						elif "0" in ccs_Freq and float(ccs_Freq) == 0:
							count0CCS += 1
							
					if Po in NGS_P_S_FreqD[CCS_Samp]:
						ngs_Freq = NGS_P_S_FreqD[CCS_Samp][Po]
						if ngs_Freq == "NA":
							countNA += 1
						if "0" in ngs_Freq and float(ngs_Freq) == 0:
							count0NGS += 1
							#print(ngs_Freq)
							#print("NGS")

					if countNA == 0 and count0NGS + count0CCS != 2:
						Os.append("\t".join([Po,ccs_Freq,ngs_Freq]))
						#print(count0NGS, count0CCS)
						#print("\t".join([Po,ccs_Freq,ngs_Freq]))
				O = "\n".join(Os)
				print(O)
				out_FO.write(O)
	out_FO.close()







if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-C","--CCS_table_P",
					  dest = "CCS_table_P",
					  default = "",
					  metavar = "path",
					  help = "path of CCS .  [required]")

	parser.add_option("-N","--NGS_table_P",
					  dest = "NGS_table_P",
					  default = "",
					  metavar = "path",
					  help = "path of NGS .  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output path to save the Freq list of NGS and CCS seq.  [required]")



	(options,args) = parser.parse_args()
	CCS_table_P   = os.path.abspath(options.CCS_table_P)
	NGS_table_P	 = os.path.abspath(options.NGS_table_P)
	out_P		 = os.path.abspath(options.out_P)



	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()





