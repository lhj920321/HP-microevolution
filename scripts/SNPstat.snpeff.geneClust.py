
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 


def Read_pan_genes(pan_genes_csv,pangenesLen):
	pangenomels = open(pan_genes_csv).readlines()
	HeaderL = pangenomels[0].split("\n")[0].split("\",\"")
	lie = 0
	lieD = {}
	for Header in HeaderL:
		if "P9-E-t" in Header:
			Header= Header.split("\"")[0]
		lieD[lie] = Header
		lie += 1


	SampGenePref = {}
	for pangenomel in pangenomels:
		if "Gene" not in pangenomel and pangenomel != "\n" :
			lL = pangenomel.split("\n")[0].split("\",\"")
			#geneName = lL[0]
			geneID = lL[0].split("\"")[1]
			for INdex in range(14,len(lL)):
				sampgeneID = lieD[INdex]
				Gid = lL[INdex]
				if Gid != "":
					if len(Gid.split("\t")) != 1 :
						for GID in Gid.split("\t"):
							SampGenePref[GID] = geneID
					else:
						SampGenePref[Gid] = geneID
	
				#else:
					#break

	return SampGenePref


'''
	out_sampF = out_P + "/Samps_SNPvcfSnpeff.stat.txt"
	if ( os.path.exists(out_sampF)):
		os.remove(out_sampF)
	out_P_O=open(out_sampF,'a')
	#out_P_O.write(LenStat + "\n")
	out_P_O.close()	

'''


def iSNV_SNP_tableRead(iSNV_SNP_table):
	iSNVSNPtablels = open(iSNV_SNP_table,'r').readlines()
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
	return lieD,lsD



def samp_Freq_Dic(samples,lieD,lsD):
	samp_FreqDic = {}
	sampNA_PosiDic = {}
	for sample in samples:
		Lie = lieD[sample]
		samp_FreqDic[sample] = {}
		sampNA_PosiDic[sample] = []
		for Posi in lsD:
			Freq = lsD[Posi][Lie]
			if Freq != "NA" :
				if Freq != "NO":				
					Freq = float(lsD[Posi][Lie])
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



			


def pan_genome_genes(pan_genesF):
	pangenesLen = {}
	pan_genesFls = open(pan_genesF).readlines()
	for pan_genesFl in pan_genesFls:
		if ">" in pan_genesFl:
			geneID = pan_genesFl.split("\n")[0].split(" ")[1]
			seq = ''
		else:
			seq += pan_genesFl.split("\n")[0]
		pangenesLen[geneID] = len(seq)
	#print(pangenesLen.keys())
	return pangenesLen



def readVcf(VcfP,persNAPosiD,SampGenePref):
	files = os.listdir(VcfP)
	GeneSNPCountDic = {}
	for File in files:
		VcfF = VcfP + "/" + File
		PersID = File.split("-")[0]
		PersNAPosi = persNAPosiD[PersID]
		SampID = File.split(".minFreq0.02.snpEffAnno.vcf")[0].split("_2_")[0]		
		refID = File.split(".minFreq0.02.snpEffAnno.vcf")[0].split("_2_")[1]
		#Geneid = SampGenePref[refID]
		
		GeneSNPCountDic[SampID] = {}

		VcfFls = open(VcfF).readlines()
		for VcfFl in VcfFls:
			if "#" not in VcfFl:
				Tags = VcfFl.split("\n")[0].split("\t")
				Cite = Tags[1]
				Info = Tags[7]
				Freq = float(Tags[9].split(":")[5].split("%")[0])
				if Cite not in PersNAPosi and Freq >= minFreqofSNP * 100:
					InfoTags = Info.split(";")
					for InfoTag in InfoTags:
						if "ANN=" in InfoTag:
							Tgs = InfoTag.split("|")
							Gene = Tgs[3]
							if Gene in SampGenePref:
								GeneID = SampGenePref[Gene]
							else:
								GeneID = Gene

							Eff = Tgs[1]
							if GeneID not in  GeneSNPCountDic[SampID]:
								GeneSNPCountDic[SampID][GeneID] = {}
							if Eff not in GeneSNPCountDic[SampID][GeneID]:
								GeneSNPCountDic[SampID][GeneID][Eff] = 0
							GeneSNPCountDic[SampID][GeneID][Eff] += 1

	#print(GeneSNPCountDic)

	return GeneSNPCountDic					


def out(GeneSNPCountDic,pangenesLen):
	SNP_stat_out = []
	clustGenes = {}
	allclustGenes = {}
	for Samp in GeneSNPCountDic:
		Samp_count = {}
		if Samp not in clustGenes:
			clustGenes[Samp] = {}
		for Gene in GeneSNPCountDic[Samp]:
			GeneSNPnum = 0
			for Eff in GeneSNPCountDic[Samp][Gene]:
				if Eff not in Samp_count:
					Samp_count[Eff] = 0
				Samp_count[Eff] += GeneSNPCountDic[Samp][Gene][Eff]
				GeneSNPnum += GeneSNPCountDic[Samp][Gene][Eff]
				#print(GeneSNPnum)
			if GeneSNPnum >= 3:
			#if "missense_variant" in GeneSNPCountDic[Samp][Gene] and GeneSNPCountDic[Samp][Gene]["missense_variant"] >= 3:
				clustGenes[Samp][Gene] = GeneSNPnum
				#print("hhshs")
				if Gene not in allclustGenes:
					allclustGenes[Gene] = ''

		synNum,MissNum,nonCodingNum,StopNum = 0,0,0,0
		if "synonymous_variant" in Samp_count:
			synNum = Samp_count["synonymous_variant"]
		
		if "missense_variant" in Samp_count:
			MissNum = Samp_count["missense_variant"]

		if "upstream_gene_variant" in Samp_count :
			nonCodingNum += Samp_count["upstream_gene_variant"]
		if "downstream_gene_variant" in Samp_count :
			nonCodingNum += Samp_count["downstream_gene_variant"]
		
		for efftype in Samp_count:
			if "stop" in efftype:
				StopNum = Samp_count[efftype]


		SNPCount = "\t".join([Samp,str(synNum),str(MissNum),str(nonCodingNum),str(StopNum)])
		SNP_stat_out.append(SNPCount)
		print(SNPCount)
	sampLst = []
	for Samp in GeneSNPCountDic:
		sampLst.append(Samp)

	#print(allclustGenes)
	sampLst.sort()
	#print(sampLst)
	outLs,outDensity = [],[]
	for clustGene in allclustGenes:
		if "group_" not in clustGene:
			outL = [clustGene]
			outD =  [clustGene]
			for Samp in sampLst:
				if clustGene in clustGenes[Samp]:
					clustNum = str(clustGenes[Samp][clustGene])
					clustDensity = str(round(clustGenes[Samp][clustGene]/pangenesLen[clustGene],4))
				else:
					clustNum = "-"
					clustDensity = "-"
				outL.append(clustNum)
				outD.append(clustDensity)
			outl = "\t".join(outL)
			outDenl = "\t".join(outD)
			outLs.append(outl)
			outDensity.append(outDenl)
	ls = "\n".join(outLs)
	Dls = "\n".join(outDensity)
	#print(ls)


	SNP_stat = "\n".join(SNP_stat_out)
	StatoutF = out_P + "/SNP.min"+ str(minFreqofSNP) + "_ann.stat.txt"
	if (os.path.exists(StatoutF)):
		os.remove(StatoutF)
	out_P_O=open(StatoutF,'a')
	out_P_O.write("\t".join(["samp","synNum","MissNum","nonCodingNum","stopNum"]) + "\n")
	out_P_O.write(SNP_stat + "\n")
	out_P_O.close()	
			



	ClustoutF = out_P + "/SNP.min"+ str(minFreqofSNP) + "_ann.clustGenes.stat.txt"
	if (os.path.exists(ClustoutF)):
		os.remove(ClustoutF)
	out_P_O=open(ClustoutF,'a')
	out_P_O.write("\t".join(["gene"] + sampLst) + "\n")
	out_P_O.write(ls + "\n")
	out_P_O.close()	
			


	ClustDensoutF = out_P + "/SNP.min"+ str(minFreqofSNP) + "_ann.clustGenes.density.stat.txt"
	if (os.path.exists(ClustDensoutF)):
		os.remove(ClustDensoutF)
	out_P_O=open(ClustDensoutF,'a')
	out_P_O.write("\t".join(["gene"] + sampLst) + "\n")
	out_P_O.write(Dls + "\n")
	out_P_O.close()	

	


def main():
	pan_genesF = roaryP + "/pan_genome_reference.fa"
	pan_genes_csv = roaryP + "/gene_presence_absence.csv"
	pangenesLen=pan_genome_genes(pan_genesF)
	SampGenePref = Read_pan_genes(pan_genes_csv,pangenesLen)

	persNAPosiD = {}
	
	for PersonIdx in range(1,26):
		if PersonIdx != 13 and PersonIdx != 14 :
			person = "P" + str(PersonIdx)
			iSNV_SNP_table = iSNV_tableP + "/" + person + "/all.iSNV_with_SNP.pyResults.txt"
			lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
			sampIDs = []
			for sampID in list(lieD.keys()):
				if "#Posi" not in sampID and "snv/SNP" not in sampID and sampID not in ["P24-C_sC-j_2_P24-C_sC-t","P8-E-t_2_P8-E-x"]:
					sampIDs.append(sampID)
			samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(sampIDs,lieD,lsD)
			personNAPosi = Person_NA_Posi(sampIDs,sampNA_PosiDic)
			persNAPosiD[person] = personNAPosi


	GeneSNPCountDic = readVcf(vcfP,persNAPosiD,SampGenePref)
	out(GeneSNPCountDic,pangenesLen)

	#out_P_O=open(out_P,'a')
	#out_P_O.write(stat + "\n")
	#out_P_O.close()	







if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--roaryP",
					  dest = "roaryP",
					  default = "",
					  metavar = "path",
					  help = "Path of roary output files .  [required]")

	parser.add_option("-v","--vcfP",
					  dest = "vcfP",
					  default = "",
					  metavar = "path",
					  help = "Path of snpeff vcf files .  [required]")

	parser.add_option("-s","--iSNV_tableP",
					  dest = "iSNV_tableP",
					  default = "",
					  metavar = "path",
					  help = "Path of iSNV table .  [required]")
	parser.add_option("-m","--minFreqofSNP",
					  dest = "minFreqofSNP",
					  default = "0.95",
					  metavar = "float",
					  help = "min Freq of SNP .  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "file",
					  help = "output stat file Path of snv Freq distribution.  [required]")



	(options,args) = parser.parse_args()
	roaryP  = os.path.abspath(options.roaryP)
	vcfP    = os.path.abspath(options.vcfP)
	iSNV_tableP    = os.path.abspath(options.iSNV_tableP)
	out_P	= os.path.abspath(options.out_P)
	minFreqofSNP	= float(options.minFreqofSNP)


	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()





