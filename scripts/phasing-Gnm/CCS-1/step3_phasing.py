#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os,time,re 
from optparse import OptionParser
import numpy as np


def readMergeReadSTypes(MergeReadSTypesF):
	SeqD,SeqDep = {},{}
	for MergeReadSTypesFl in open(MergeReadSTypesF).readlines():
		Tags = MergeReadSTypesFl.split("\n")[0].split("\t")
		supNum = Tags[0]
		PosiS = Tags[1].split("-")[1]
		PosiE = Tags[1].split("-")[2]
		Seq = ''
		for PoBs in  Tags[2].split("_")[:-1]:
			Po = re.findall("\d+",PoBs)[0]
			#print(Po)
			Bs = PoBs[-1]
			Seq += Bs
			#print(Bs)
		SeqD[PosiS + "_" + PosiE + "_" + Seq] = Seq
		SeqDep[PosiS + "_" + PosiE + "_" + Seq] = supNum
		#print(PosiS)
		#print(Seq)
	return SeqD,SeqDep


def readPosi(filt_posi_vcfF):	
	filtPosiLst,filtPosiD = [],{}
	filt_posi_ls = open(filt_posi_vcfF,'r').readlines()
	Count = 0
	for filt_posi_l in filt_posi_ls:
		filtPosiLst.append(int(filt_posi_l.split("\t")[1]))
		filtPosiD[int(filt_posi_l.split("\t")[1])] = Count
		Count += 1
	return filtPosiLst,filtPosiD



def phasing_search(filtPosiLst,filtPosiD,SeqD,SeqDep,startLociNum):
	StarSidx = 0
	#for searchIndxS in range(len(filtPosiLst)-10):
	RightSeq,outPSeq = {},{}
	for searchIndxS in range(20):
	#for searchIndxS in range(len(filtPosiLst)-startLociNum):
		if searchIndxS  >= StarSidx:
			print("search start index :" + str(searchIndxS))
			#Ct = 0
			TreeSeqNumD,PhasingFlagD = {},{}
			for searchIndxE in range(searchIndxS + startLociNum,20):
			#for searchIndxE in range(searchIndxS + startLociNum,len(filtPosiLst)):
				#Ct += 1
				TreeSeqNumD[searchIndxE] = 0
				print("search end index :" + str(searchIndxE-1))
				searchPosS = filtPosiLst[searchIndxS]
				searchPosE = filtPosiLst[searchIndxE]
				SearchLen  = searchIndxE - searchIndxS 
				searchPoss = filtPosiLst[searchIndxS:searchIndxE]
				#print(searchIndxS)
				#print(searchIndxE)		
				#print("search Posi :")
				#print(searchPoss)

				TreeSeqNum = 0
				seqBs = {}

				for SeqSE in SeqD:
					RS = int(SeqSE.split("_")[0])
					RE = int(SeqSE.split("_")[1])
					if RS <= searchPoss[0] and RE >= searchPoss[0]:
						RS = searchPoss[0]
					if RE >= searchPoss[-1] and RS <= searchPoss[-1]:
						RE = searchPoss[-1]
						
					RSidx = filtPosiD[RS]
					REidx = filtPosiD[RE]
					SEQ = SeqD[SeqSE]
					Len = len(SeqD[SeqSE])	
					allPosis= filtPosiLst[RSidx:REidx+1]
					#print("shishishishishishihsi")
					if REidx >= searchIndxE:
						CovLen = searchIndxE - RSidx
					else:
						CovLen = REidx - RSidx 

					if RS > searchPosE or RE < searchPosS:
						continue
					elif CovLen >= minCov * SearchLen:
							#print("ososoosososoo")
							#print(CovLen)
							if searchIndxS not in RightSeq:
								RightSeq[searchIndxS] = {}
							if searchIndxE not in RightSeq[searchIndxS]:
								RightSeq[searchIndxS][searchIndxE] = {}
							RightSeq[searchIndxS][searchIndxE][SeqSE] = SeqDep[SeqSE]
							#print(RightSeq)
							#print(SeqD[SeqSE])
							for allPoIndex in range(len(allPosis)):
								allPo = allPosis[allPoIndex]
								Bs = SeqD[SeqSE][allPoIndex]
								if Bs != "-" :
									if allPo not in seqBs:
										seqBs[allPo] = 0
									seqBs[allPo] += 1

						#if CovLen >= minCov * SearchLen:
							TreeSeqNumD[searchIndxE] += 1



				#print(len(searchPoss))
				#print(seqBs)
				if len(seqBs) == len(searchPoss):
					Flag = "Can"
				else:
					Flag = "Cannot"

				PhasingFlagD[searchIndxE] = Flag
				if searchIndxE >= searchIndxS + startLociNum + 1:
					#print("------------")
					#print(PhasingFlagD[searchIndxE])
					#print(PhasingFlagD[searchIndxE-1])
					#print("oooooooooooo")
					if PhasingFlagD[searchIndxE-1] == "Can" and PhasingFlagD[searchIndxE] == "Cannot": 
					#if PhasingFlagD[searchIndxE] == "Cannot": 
						StarSidx = searchIndxE -1
						if searchIndxS not in outPSeq:
							#print("osososoosoosososososoososososoososo")
							outPSeq[str(searchIndxS) + "_" + str(filtPosiLst[searchIndxS]) + "_" + str(filtPosiLst[searchIndxE-1])] = RightSeq[searchIndxS][searchIndxE-1]
							#print(outPSeq)
						
						print("break !!!")
						#print(RightSeq[searchIndxS][searchIndxE-1])
						break
					if PhasingFlagD[searchIndxE-1] == "Cannot" and PhasingFlagD[searchIndxE] == "Cannot":
						break
				#if searchIndxS == 30 and PhasingFlagD[searchIndxE] == "Can":
					##if searchIndxE == len(filtPosiLst):
					#if searchIndxS not in outPSeq:
						#outPSeq[str(searchIndxS) + "_" + str(filtPosiLst[searchIndxS])] = RightSeq[searchIndxS][searchIndxE-1]
						#print(outPSeq)
						

				
	#print(outPSeq)
	Lenslst,containPosiNum,Outsummaryls = [],0,[]
	for RegionID in outPSeq:
		#endPo = list(outPSeq[RegionID].keys())[0].split("_")[1]
		startPo = RegionID.split("_")[1]
		endPo = RegionID.split("_")[2]
		Sidx = filtPosiD[int(startPo)]
		Eidx = filtPosiD[int(endPo)]
		PosiNum = Eidx - Sidx 
		RgAllPosis = filtPosiLst[Sidx:Eidx]
		#print(RgAllPosis)
		containPosiNum += PosiNum
		Lenslst.append(int(endPo)-int(startPo))
		outF = outP + "/S_" + startPo + "-E_" + endPo  + "-L." + str(int(endPo)-int(startPo)) + "-num_" + str(PosiNum) + ".fasta"
		Outsummaryls.append("\t".join([startPo,endPo,str(int(endPo)-int(startPo)),str(PosiNum)]))
		Ols = []
		outID_Dic = {}
		for RegionSeqID in outPSeq[RegionID]:
			rgRdS = RegionSeqID.split("_")[0]
			rgRdE = RegionSeqID.split("_")[1]
			RdSeq = RegionSeqID.split("_")[2]
			rgRdSidx = filtPosiD[int(rgRdS)]
			rgRdEidx = filtPosiD[int(rgRdE)]
			rgRdPosiNum = rgRdEidx - rgRdSidx +1
			Depth = SeqDep[RegionSeqID]
			RgRdAllPosis = filtPosiLst[rgRdSidx:rgRdEidx+1]
			RgRdAllPoBases = {}
			for RgRdAllPoIdx in range(len(RgRdAllPosis)):
				base = RdSeq[RgRdAllPoIdx]
				RgRdAllPo = RgRdAllPosis[RgRdAllPoIdx]
				RgRdAllPoBases[RgRdAllPo] = base

			#print(RgRdAllPoBases)
			outseq = ''
			for RgAllPosi in RgAllPosis:
				#print(RgAllPosi)
				if RgAllPosi in RgRdAllPosis:
					Base = RgRdAllPoBases[RgAllPosi]
				else:
					Base = '-'
				outseq += Base
			outID = ">S." + rgRdS + "-E." + rgRdE + "-L." + str(int(rgRdE)-int(rgRdS))  + ".iSNVPosiNum_" + str(rgRdPosiNum) + ".D_" + Depth
			if outID not in outID_Dic:
				outID_Dic[outID] = 0
			else:
				outID_Dic[outID] += 1
			Ols.append(outID + "--" + str(outID_Dic[outID]) )  # + str(SeqDep[])
			Ols.append(outseq)
		O = "\n".join(Ols)



		if (os.path.exists(outF)) :
			os.remove(outF)
		ONT_FO = open(outF,'a')
		ONT_FO.write(O+"\n")
		ONT_FO.close()


##summary
		StatF = outP + "/Phasing_summary.txt"
		Len_mean = round(np.mean(Lenslst),2)
		Len_sum = np.sum(Lenslst)
		Lenslst.sort()
		Len_median = np.median(Lenslst)
		summary="\t".join([str(len(Lenslst)),str(Len_mean),str(Len_median),str(Len_sum),str(containPosiNum),str(len(filtPosiLst)),str(float(round(containPosiNum *100/len(filtPosiLst),2)))])
		Head = "\t".join(["contigNum","MeanLen","MedianLen","Len_sum","PhasedLociNum","allReliaLociNum","percentLociNum(%)"])
		if (os.path.exists(StatF)) :
			os.remove(StatF)
		ONT_FO = open(StatF,'a')
		ONT_FO.write(Head+"\n")
		ONT_FO.write(summary+"\n")
		ONT_FO.write("\n")
		ONT_FO.write("\n")

		ONT_FO.write("Phased Region summary : "+ "\n")
		ONT_FO.write("\n")
		
		RegionSummaryH="\t".join(["StartLoci","EndLoci","RegionLen","PhasedLociNum"])
		RegionSummaryls = "\n".join(Outsummaryls)

		ONT_FO.write(RegionSummaryH+"\n")
		ONT_FO.write(RegionSummaryls+"\n")

		ONT_FO.close()






def main():
	startTime = time.time()
	################################################################################################################################
	# Check input files
	################################################################################################################################

	filtPosiLst,filtPosiD = readPosi(filt_posi_vcfF)
	print("all Posi Num  :  " + str(len(filtPosiLst)))
	print(MergeReadSTypesF)
	SeqD,SeqDep = readMergeReadSTypes(MergeReadSTypesF)
	#print(SeqD)


	phasing_search(filtPosiLst,filtPosiD,SeqD,SeqDep,startLociNum)



##log 
	if (os.path.exists(out_line_file)) :
		os.remove(out_line_file)
	ONT_FO = open(out_line_file,'a')
	ONT_FO.write("step:3"+"\n")
	ONT_FO.close()


	endTime = time.time()
	sys.stdout.write("Total time taken: "+str(endTime-startTime)+" seconds\n")








if __name__ == "__main__":
	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)


	#parser.add_option("-c","--config_F",
					  #dest = "config_F",
					  #default = "",
					  #metavar = "file",
					  #help = "Config file [required]")
	parser.add_option("-f","--filt_posi_vcfF",
					  dest = "filt_posi_vcfF",
					  default = "",
					  metavar = "string",
					  help = "filt posi File [required]")

	parser.add_option("-F","--MergeReadSTypesF",
					  dest = "MergeReadSTypesF",
					  default = "",
					  metavar = "file",
					  help = "out path  [required]")

	parser.add_option("-o","--outP",
					  dest = "outP",
					  default = "",
					  metavar = "path",
					  help = "out path  [required]")

	parser.add_option("-l","--out_line_file",
					  dest = "out_line_file",
					  default = "",
					  metavar = "file",
					  help = "out_line_file . [required]")
	parser.add_option("-s","--startLociNum",
					  dest = "startLociNum",
					  default = "",
					  metavar = "int",
					  help = "start search Loci Num . [required]")
	parser.add_option("-m","--minLociNum_forOut",
					  dest = "minLociNum_forOut",
					  default = "",
					  metavar = "file",
					  help = "min Loci Num for Output . [required]")
	parser.add_option("-c","--minCov",
					  dest = "minCov",
					  default = "0.9",
					  metavar = "float",
					  help = "min Coverage of the read in the phasing region  . [required]")


	(options,args) = parser.parse_args()

	#config_F	          = os.path.abspath(options.config_F)
	filt_posi_vcfF		  = os.path.abspath(options.filt_posi_vcfF)
	MergeReadSTypesF      = os.path.abspath(options.MergeReadSTypesF)	
	outP                  = os.path.abspath(options.outP)	
	out_line_file	      = os.path.abspath(options.out_line_file)
	startLociNum	      = int(options.startLociNum)
	minLociNum_forOut	  = options.minLociNum_forOut
	minCov                = float(options.minCov)

	main()













