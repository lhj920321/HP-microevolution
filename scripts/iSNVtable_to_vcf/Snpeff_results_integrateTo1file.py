#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os,time 
from optparse import OptionParser

#This script is use to organize all the samples' mutations annotions by the SnpEff program into a single file.
################################################################################################################################
#
# For help as a standalone program type: python Snpeff_results_integrateTo1file -h
#
# Examples:
#	python Snpeff_results_integrateTo1file -i  $outsnpeffP   -o  ./allSnpEff.txt
#	#parameters:
#		-i  input path,output path of iSNV_calling.sh program."iSNV_with_SNP.all.txt","iSNV_info.txt"and"SNP_info.txt" are required in the input path.
#		-o  output file, Containing annotated information for all sample mutations using snpeff software.
#
#	
################################################################################################################################

def charRank_geneRegion(gene,posi):
	refGeneDic = {"orf1a":"266..13468","orf1b":"13468..21555","S":"21563..25384"\
	,"ORF3a":"25393..26220","E":"26245..26472","M":"26523..27191",\
	"ORF6":"27202..27387","ORF7a":"27394..27759","ORF8":"27894..28259",\
	"N":"28274..29533","ORF10":"29558..29674",}
	for refGene in refGeneDic:
		geneIDStrt = int(refGeneDic[refGene].split("..")[0])
		geneIDEnd = int(refGeneDic[refGene].split("..")[1])
		if int(posi)  >= geneIDStrt  and int(posi)  <= geneIDEnd:
			charRank = (int(posi)-geneIDStrt )%3 + 1
			break

	return charRank


def sampSnpeffStat(sample,sampSnpeffF,outAllSampF):
	outAllSampFO = open(outAllSampF,'a')
	sampSnpeffFls = open(sampSnpeffF).readlines()
	for sampSnpeffFl in sampSnpeffFls:
		if "#" not in sampSnpeffFl and sampSnpeffFl != "\n":
			outLst = []
			linekey = sampSnpeffFl.split("\n")[0].split("\t")
			posi = 	linekey[1]		
			ref = 	linekey[3]
			alle = 	linekey[4]
			ann = 	linekey[7]
			annkeys = ann.split(";")[5].split("|")
			Annotation = annkeys[1]
			Gene_Name = annkeys[3]
			Feature_Type = annkeys[5]
			Transcript_BioType = annkeys[7]
			
			if Annotation not in ["upstream_gene_variant",'downstream_gene_variant']:
				rank = annkeys[8].split("/")[0]
				rankall = annkeys[8].split("/")[1]
				cDNA_pos = annkeys[11].split("/")[0]
				cDNA_posall = annkeys[11].split("/")[1]
				AA_pos = annkeys[13].split("/")[0]
				AA_posall = annkeys[13].split("/")[1]
				charRank = charRank_geneRegion(Gene_Name,posi)
			else:
				rank = ''
				rankall = ''
				cDNA_pos = ''
				cDNA_posall = ''
				AA_pos = ''
				AA_posall = ''
				charRank = ''


			HGVS_c = annkeys[9]
			HGVS_p = annkeys[10]

			#CDS_pos = annkeys[12]

			Freq = linekey[-1].split(":")[-1].split("%")[0]
			
			outLst.append(sample)
			outLst.append(posi)
			outLst.append(ref)
			outLst.append(alle)
			outLst.append(Annotation)
			outLst.append(Gene_Name)
			outLst.append(Feature_Type)
			outLst.append(Transcript_BioType)
			outLst.append(rank)
			outLst.append(rankall)
			outLst.append(HGVS_c)
			outLst.append(HGVS_p)
			outLst.append(cDNA_pos)
			outLst.append(cDNA_posall)
			outLst.append(AA_pos)
			outLst.append(AA_posall)
			outLst.append(str(charRank))
			outLst.append(Freq)
			out  = "\t".join(outLst)
			outAllSampFO.write(out + "\n")
	outAllSampFO.close()


def header():
	outLst = []
	outLst.append('sample')
	outLst.append('posi')
	outLst.append('ref')
	outLst.append('alle')
	outLst.append('Annotation')
	outLst.append('Gene_Name')
	outLst.append('Feature_Type')
	outLst.append('Transcript_BioType')
	outLst.append('rank')
	outLst.append('rankall')
	outLst.append('HGVS_c')
	outLst.append('HGVS_p')
	outLst.append('cDNA_pos')
	outLst.append('cDNA_posall')
	outLst.append('AA_pos')
	outLst.append('AA_posall')
	outLst.append('charRank')
	outLst.append('Freq%')
	outHead = "\t".join(outLst)
	return outHead






def main():	
	startTime = time.time()

	outHead = header()
	outAllSampFO = open(outAllSampF,'a')
	outAllSampFO.write(outHead + "\n")
	outAllSampFO.close()

	filesLst = os.listdir(snpeffP)
	count = 1
	for SnpeffF in filesLst:
		sampSnpeffF = snpeffP + "/" + SnpeffF
		print "file no." + str(count) + ":   " + SnpeffF
		sample = sampSnpeffF.split("/")[-1].split(".iSNV")[0]
		print(sample)
		sampSnpeffStat(sample,sampSnpeffF,outAllSampF)
		count +=1

	endTime = time.time()
	print 
	print "Total samples: " + str(count -1)
	sys.stdout.write("Total time taken: "+str(endTime-startTime)+" seconds\n")


if __name__ == "__main__":
	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--inputpath",
	                  dest = "inputpath",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to input vcf files annotioned by SnpEff program [required]")

	parser.add_option("-o","--outputfile",
	                  dest = "outputfile",
	                  default = "",
	                  metavar = "file",
	                  help = "output file, Containing annotated information for all sample mutations using snpeff software. [required]")

	(options,args) = parser.parse_args()


	snpeffP      = os.path.abspath(options.inputpath)
	outAllSampF  = os.path.abspath(options.outputfile)


	if (os.path.exists(outAllSampF)):
		os.remove(outAllSampF)


	#FREQ_threshold = 1
	#baseCovThreshold = 10

	

	main()












