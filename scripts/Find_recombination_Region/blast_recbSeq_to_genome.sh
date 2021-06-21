
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
	PersonDic = {}
	lieD = {}

	for Header in HeaderL:
		if "P9-E-t" in Header:
			Header= Header.split("\"")[0]

		lieD[lie] = Header
		lie += 1
		if lie > 14 and  Header != "P11-E-j.ONT" and Header != "P8-E-t" and Header != "P24-C_sC-j":			 
			Person = Header.split("-")[0]
			if Person not in PersonDic:
				PersonDic[Person] = {}
			PersonDic[Person][Header] = ''
	print(PersonDic)

	Samp_pan_genesD = {}
	Person_panGene_count = {}
	linecout = 0
	allPerson_panGeneD = {}
	
	for pangenomel in pangenomels:
		linecout += 1
		if "Gene" not in pangenomel and pangenomel != "\n" :
			lL = pangenomel.split("\n")[0].split("\",\"")
			geneID = lL[0].split("\"")[1]
			for INdex in range(14,len(lL)):
				sampID = lieD[INdex]
				if sampID != "P11-E-j.ONT" and sampID != "P8-E-t" and sampID != "P24-C_sC-j":			 
					persID = sampID.split("-")[0]
					Gid = lL[INdex]
					if sampID == "P9-E-t":
						Gid = lL[INdex].split("\"")[0]
					if Gid != '':
						if sampID not in Samp_pan_genesD:
							Samp_pan_genesD[sampID] = {}
						Samp_pan_genesD[sampID][geneID] = ''
						if persID not in Person_panGene_count:
							Person_panGene_count[persID] = {}					
						if geneID not in Person_panGene_count[persID]:
							Person_panGene_count[persID][geneID] = []
						Person_panGene_count[persID][geneID].append(sampID)
						if geneID not in allPerson_panGeneD:
							allPerson_panGeneD[geneID] = ''
	allPerson_panGeneLen = 0
	for GE in allPerson_panGeneD:
		allPerson_panGeneLen += pangenesLen[GE]					

	sampUniqGene = {}
	sampUniqGeneLen = {}
	outs = []
	outs.append("\t".join(["person","diffGeneNum","diffGene_gene_Len"]))
	for person in PersonDic:
		person_diff_geneC = 0
		diffGene_gene_Len = 0

		for GENE in  Person_panGene_count[person]:
			if len(Person_panGene_count[person][GENE]) != 0 and len(Person_panGene_count[person][GENE]) != len(PersonDic[person]):
				person_diff_geneC += 1
				diffGene_gene_Len += pangenesLen[GENE]	
			if len(Person_panGene_count[person][GENE]) != 0 and len(Person_panGene_count[person][GENE]) == 1:
				uniqGeneSamp = Person_panGene_count[person][GENE][0]
				if uniqGeneSamp not in sampUniqGene:
					sampUniqGene[uniqGeneSamp] = 0
					sampUniqGeneLen[uniqGeneSamp] = 0
				sampUniqGene[uniqGeneSamp] += 1
				sampUniqGeneLen[uniqGeneSamp] += pangenesLen[GENE]

		outs.append("\t".join([person,str(person_diff_geneC),str(diffGene_gene_Len)]))
	outPs = "\n".join(outs)
	#print(outPs)

	LenOut = []
	LenOut.append("\t".join(["sample","pan_gene_Num","pan_gene_Len","uniqGeneNum","uniqGeneLen"]))

	print(sampUniqGene)
	for samp in Samp_pan_genesD:
		panGenLen = 0
		for geneid in Samp_pan_genesD[samp]:
			panGenLen += pangenesLen[geneid]	
		print(samp)	
		if samp in sampUniqGene:		
			LenOut.append("\t".join([samp,str(len(Samp_pan_genesD[samp])),str(panGenLen),str(sampUniqGene[samp]),str(sampUniqGeneLen[samp])]))
		else:
			LenOut.append("\t".join([samp,str(len(Samp_pan_genesD[samp])),str(panGenLen),'0','0']))

	LenStat = "\n".join(LenOut)
			
	print(outPs) 
	print(LenStat) 

	stat = " ".join(["all samples pan-genome:",str(len(allPerson_panGeneD)),str(allPerson_panGeneLen)])
	print(stat)


	out_F = out_P + "/Persons_pan-gene.stat.txt"
	if ( os.path.exists(out_F)):
		os.remove(out_F)
	out_P_O=open(out_F,'a')
	out_P_O.write(stat + "\n")
	out_P_O.write(outPs + "\n")
	out_P_O.close()	

	out_sampF = out_P + "/Samps_pan-gene.stat.txt"
	if ( os.path.exists(out_sampF)):
		os.remove(out_sampF)
	out_P_O=open(out_sampF,'a')
	out_P_O.write(LenStat + "\n")
	out_P_O.close()	




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





def main():
	pan_genesF = roaryP + "/pan_genome_reference.fa"
	pan_genes_csv = roaryP + "/gene_presence_absence.csv"
	pangenesLen=pan_genome_genes(pan_genesF)
	Read_pan_genes(pan_genes_csv,pangenesLen)

#
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
	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "file",
					  help = "output stat file Path of snv Freq distribution.  [required]")



	(options,args) = parser.parse_args()
	roaryP   = os.path.abspath(options.roaryP)
	out_P		  = os.path.abspath(options.out_P)



	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)
		#os.remove(out_P)


	main()





