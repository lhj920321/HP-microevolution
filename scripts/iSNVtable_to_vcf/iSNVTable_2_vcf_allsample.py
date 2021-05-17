#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os,time 
from optparse import OptionParser
import revert_modules 



# The role of this script is to convert iSNV table into vcf format(VCFv4.1) that can be recognized  by other mutation annotion programs
################################################################################################################################
#
# For help as a standalone program type: python iSNVTable_2_vcf_allsample.py -h
#
# Examples:
#	    python iSNVTable_2_vcf_allsample.py  -i  ../data/   -o  $isnvTableVcfP   -r   MN908947.3
#parameters:
#			-i  input path,the path to the  output files of the program named "iSNV_calling.sh" ,these files need to be included in this path: "iSNV_with_SNP.all.txt","iSNV_info.txt"and"SNP_info.txt" .
#			-o  output path, the path to save the output vcf files.
#			-r  The reference genome ID will be output to the output files,same as the ref ID used in "iSNV_info.txt" and "SNP_info.txt" files.
#
# Notes:
#    -  Everything is calculated against the reference. 
#	
#	
################################################################################################################################


def isnvTableSamp(SAMPIDs,iSNVTableP,outvcfP,FREQ_threshold,baseCovThreshold,RefChromID):
	iSNVTable = iSNVTableP + "/"+ "iSNV_with_SNP.all.txt"
	iSNVinfoF = iSNVTableP + "/"+ "iSNV_info.txt"
	SNPinfoF = iSNVTableP + "/"+ "SNP_info.txt"

	for SAMPID in SAMPIDs:
		print "samp:  " + SAMPID
		outvcfF = outvcfP + "/" + SAMPID + ".iSNV.vcf"
		Fexit = os.path.exists(outvcfF)
		if Fexit == True:
			os.remove(outvcfF)
		header = revert_modules.vcf_Header(SAMPID) 
		outvcfFO = open(outvcfF,'a')
		outvcfFO.write(header + "\n")
		outvcfFO.close()
	
		sampPosL = revert_modules.iSNVsampPos(iSNVTable,SAMPID,FREQ_threshold)		
		snvInfoD = revert_modules.samp_posi_info(iSNVinfoF,sampPosL,SAMPID,baseCovThreshold)
		SNPInfoD = revert_modules.samp_posi_info(SNPinfoF,sampPosL,SAMPID,baseCovThreshold)

		out = revert_modules.vcfFomat(sampPosL,snvInfoD,SNPInfoD,RefChromID)
		outvcfFO = open(outvcfF,'a')
		outvcfFO.write(out + "\n")
		outvcfFO.close()


def main():
	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--inputpath",
	                  dest = "inputpath",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to input iSNV tables [required]")

	parser.add_option("-o","--outputpath",
	                  dest = "outputpath",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to save output files in [required]")
	parser.add_option("-r","--RefChromID",
	                  dest = "RefChromID",
	                  default = "",
	                  metavar = "string",
	                  help = "The chrom ID of reference genome in the output vcf files [required]")

	parser.add_option("-e","--FREQ_threshold",
	                  dest = "FREQ_threshold",
	                  default = "0.05",
	                  metavar = "",
	                  help = "Minmum frequency of iSNV for output ,[default: %default][optional]")

	parser.add_option("-m","--MinReadsSupport",
	                  dest = "MinReadsSupport",
	                  default = "10",
	                  metavar = "",
	                  help = "Minmum num of reads to support the mutation ,[default: %default][optional]")


	(options,args) = parser.parse_args()

	iSNVTableP       = os.path.abspath(options.inputpath)
	outvcfP          = os.path.abspath(options.outputpath)
	FREQ_threshold   = options.FREQ_threshold
	baseCovThreshold = options.MinReadsSupport
	RefChromID       = options.RefChromID

	print outvcfP
	if (not os.path.exists(outvcfP)):
		os.mkdir(outvcfP)

	################################################################################################################################
	# program
	################################################################################################################################
	startTime = time.time()

	iSNVTable = iSNVTableP + "/" + "iSNV_with_SNP.all.txt"
	SAMPIDs = revert_modules.sample_ID(iSNVTable)
	isnvTableSamp(SAMPIDs,iSNVTableP,outvcfP,FREQ_threshold,baseCovThreshold,RefChromID)

	endTime = time.time()
	print 
	print "Total samples: " + str(len(SAMPIDs))
	sys.stdout.write("Total time taken: "+str(endTime-startTime)+" seconds\n")


if __name__ == "__main__":
	main()












