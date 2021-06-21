
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 
def PersonRespectSamp(person):
	if person == 'P1':
		ref_samp='P1-E-j'
	if person == 'P2':
		ref_samp="P2-E-t"
	if person == 'P3':
		ref_samp="P3-E-j"
	if person == 'P4':
		ref_samp="P4-E-j"
	if person == 'P5':
		ref_samp="P5-E-d"
	if person == 'P6':
		ref_samp="P6-E-t"
	if person == 'P7':
		ref_samp="P7-E-j"
	if person == 'P8':
		ref_samp="P8-E-x"
	if person == 'P9':
		ref_samp="P9-E-j"
	if person == 'P10':
		ref_samp="P10-E-j"
	if person == 'P11':
		ref_samp="P11-E-x"
	if person == 'P12':
		ref_samp="P12-E-t"
	if person == 'P13':
		ref_samp="P13-E-j"
	if person == 'P14':
		ref_samp="P14-E-j"
	if person == 'P15':
		ref_samp="P15-C_bS-j"
	if person == 'P16':
		ref_samp="P16-C_bS-j"
	if person == 'P17':
		ref_samp="P17-C_bR-j"
	if person == 'P18':
		ref_samp="P18-C_bR-x"
	if person == 'P19':
		ref_samp="P19-C_bR-j"
	if person == 'P20':
		ref_samp="P20-C_bR-t"  ##changed
	if person == 'P21':
		ref_samp="P21-C_sL-j"
	if person == 'P22':
		ref_samp="P22-C_sL-d"  ##changed
	if person == 'P23':
		ref_samp="P23-C_sL-j"
	if person == 'P24':
		ref_samp="P24-C_sC-t"  ##changed
	if person == 'P25':
		ref_samp="P25-C_sC-t"

	return ref_samp


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



def Person_NA_Posi(samples,sampNA_PosiDic):
	personNAPosi = {}
	for sample in samples:
		sampleNAPosi = sampNA_PosiDic[sample]
		for posi in sampleNAPosi:
			if posi not in personNAPosi:
				personNAPosi[posi] = ''
	return personNAPosi


def Person_commonSNP_Posi(samples,samp_FreqDic,lsD):
	personCommonPosi = {}
	
	for posi in lsD.keys():
		count = 0
		for sample in samples:
			sampleFreqDic = samp_FreqDic[sample]
			if posi in sampleFreqDic:
				FREQ = sampleFreqDic[posi]
				if FREQ >= 0.95:
					count += 1
		if count == len(samples):
			personCommonPosi[posi] = ''
			
	return personCommonPosi




def Samps_RepeatRegion(sample,RepeatRegionP):
	SampReptReg = {}
	RepeatRegionF = RepeatRegionP + "/" + sample.split("_2_")[0] + "/" + sample.split("_2_")[0] + ".RepeatRegion.txt"
	for RepeatRegionFl in open(RepeatRegionF).readlines():
		if RepeatRegionFl != "\n":
			RepaetPosi = RepeatRegionFl.split("\t")[0]
			SampReptReg[RepaetPosi] = ''
	return SampReptReg


'''
def allSamps_recombiRgs(ClonalFrameF):
	SampsRcombRgs = {}
	for ClonalFrameFl in open(ClonalFrameF).readlines():
		if "NODE_" not in ClonalFrameFl and "Node" not in ClonalFrameFl :
			Tags = ClonalFrameFl.split("\n")[0].split("\t")
			Samp = Tags[0]
			Beg = int(Tags[1])
			End = int(Tags[2])
			for Ps  in range(Beg,End + 1):
				if Samp not in SampsRcombRgs:
					SampsRcombRgs[Samp] = {}
				SampsRcombRgs[Samp][Ps] = ''
	return SampsRcombRgs
'''







def tongjiFREQ(seqTech,sample,samp_FreqDic,out_F,personNAPosi,personCommonPosi,snvSNPcount_Dic,SampReptReg,RefSampReptReg):
	print(sample)
	#print(SampsRcombRgs.keys())
	FreqDic = samp_FreqDic[sample]
	tongji_Dic = {}
	for XX in range(0,101,step):
		tongji_Dic[XX] = 0
	snvCount,SNPCount,SNVRcombNum = 0,0,0
	SampID = sample.split("_2_")[0]
	#if  SampID in SampsRcombRgs :
		#SampRcombRg = SampsRcombRgs[SampID]
	#else:
		#SampRcombRg={}
	SNVPosiD = {}
	for Posi in FreqDic:
		if Posi not in personNAPosi and Posi not in personCommonPosi and Posi not in RefSampReptReg : #and Posi not in SampRcombRg
			Freq = FreqDic[Posi]
			for XX in range(0,101,step):
				if Freq * 100 >= XX and Freq*100 < XX+step:
					tongji_Dic[XX]  += 1

			if Freq >= minFreq and Freq < 1-minFreq:
				snvCount += 1
				SNVPosiD[Posi] = Freq
			elif Freq >= 1-minFreq :
				SNPCount += 1

		if Posi not in personNAPosi and Posi not in personCommonPosi and  (Posi in RefSampReptReg ):  #or Posi in SampRcombRg 
			SNVRcombNum += 1
	freq_flagLst=[]
	snvCount_list=[]
	outlst = []
	for XX in range(0,101,step):
		outlst.append(sample +"\t" + str(XX) + "\t" + str(tongji_Dic[XX]))
		freq_flagLst.append(XX)
		snvCount_list.append(tongji_Dic[XX])

	outlines = "\n".join(outlst) 
	out_F_O=open(out_F,'a')
	out_F_O.write(outlines)
	out_F_O.close()	
	
	out_plotF = out_F.split(".txt")[0] + ".pdf"
	plt.figure(dpi=300,figsize=(8,3))
	plt.bar(range(len(snvCount_list)), snvCount_list,color='red',tick_label=freq_flagLst)  


	plt.title(sample.split(".sort")[0])
	plt.xlabel('Freq', fontsize = 13)
	plt.ylabel('count', fontsize = 13)
	plt.xticks(rotation=90)
	plt.savefig(out_plotF)
	plt.savefig(out_plotF.split(".pdf")[0] + ".png")



## snv and SNP count 
	#out_snvSNPCountF = out_F.split(".txt")[0] + ".snvSNPcount.txt"
	#out_snvSNPCountF_O=open(out_snvSNPCountF,'a')
	#out_snvSNPCountF_O.write("\t".join(["sample","snvCount","SNPCount"]) + "\n")
	#out_snvSNPCountF_O.write("\t".join([sample,str(snvCount),str(SNPCount)]))
	#out_snvSNPCountF_O.close()	
	
	print("\t".join([sample,str(snvCount),str(SNPCount)]))
	SampNo = "{:02}".format(int(sample.split("_2_")[0].split("-")[0].split("P")[1]) )

	Sampiid = "-".join(["P" +str(SampNo), "-".join(sample.split("_2_")[0].split("-")[1:]) ] ) 
	snvSNPcount_Dic[sample] = "\t".join([seqTech,sample.split("_2_")[0].split("-")[1].split("_")[0],Sampiid,sample.split("_2_")[1],str(snvCount),str(SNPCount),str(len(personCommonPosi)),str(len(personNAPosi)),str(len(RefSampReptReg)),str(SNVRcombNum),str(len(SampReptReg))])


	##out SNV Posi and Freq
	SNVPosiOs = ["\t".join(["Sample","SeqTech","SNVPosi","Freq"])]
	for SNVPosi in SNVPosiD:
		SNVPosiOl = "\t".join([Sampiid,seqTech,SNVPosi,str(SNVPosiD[SNVPosi])])
		SNVPosiOs.append(SNVPosiOl)		
	SNVPosiO = "\n".join(SNVPosiOs)
	
	outF = out_P + "/" + sample + ".SNV.PosiAndFreq.txt"
	if (os.path.exists(outF)):
		os.remove(outF)
	outF_O=open(outF,'a')
	outF_O.write(SNVPosiO)
	outF_O.close()	


	return snvSNPcount_Dic






def main():
	lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
	iSNVtable_lst = list(lieD.keys())
	#SampsRcombRgs = allSamps_recombiRgs(ClonalFrameF)
	sampIDs = []
	for sampID in iSNVtable_lst:
		if "#Posi" not in sampID and "snv/SNP" not in sampID : #and "P11-E-t" not in sampID
			sampIDs.append(sampID)


	samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(sampIDs,lieD,lsD)
	personNAPosi = Person_NA_Posi(sampIDs,sampNA_PosiDic)
	personCommonPosi = Person_commonSNP_Posi(sampIDs,samp_FreqDic,lsD)
	snvSNPcount_Dic = {}

	for sample in sampIDs:
		person = sample.split("-")[0]
		ref_Samp = PersonRespectSamp(person)
		SampReptReg = Samps_RepeatRegion(sample,RepeatRegionP)
		RefSampReptReg = Samps_RepeatRegion(ref_Samp,RepeatRegionP)
		out_F = out_P + "/" + sample + ".iSNVpy.Freq_distribu.txt"
		snvSNPcount_Dic = tongjiFREQ(seqTech,sample,samp_FreqDic,out_F,personNAPosi,personCommonPosi,snvSNPcount_Dic,SampReptReg,RefSampReptReg)

	person = iSNV_SNP_table.split("/")[-2]
	out_snvSNPCountF = out_P + "/" + person + ".snvSNPcount.txt"
	out_snvSNPCountF_O=open(out_snvSNPCountF,'a')
	out_snvSNPCountF_O.write("\t".join(["seqTech","group","sample","RespSamp","Filted_snvCount(" + str(minFreq) + "-" + str(1-minFreq) + ")" ,"Filted_SNPCount(>=" + str(1-minFreq) + ")","personCommonCitesNum","personNAPosi","RefSampRepeatPosiNum","SNVinRefRepeat","SampRepeatPosiNum",]) + "\n")
	for sampleID in snvSNPcount_Dic:	
		out_snvSNPCountF_O.write(snvSNPcount_Dic[sampleID] + "\n")
	out_snvSNPCountF_O.close()	





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

	parser.add_option("-R","--RepeatRegionP",
					  dest = "RepeatRegionP",
					  default = "",
					  metavar = "path",
					  help = "path of genome Repeat Region file .  [required]")

	parser.add_option("-m","--minFreq",
					  dest = "minFreq",
					  default = "0",
					  metavar = "float",
					  help = "min Freq of snv to calculate (0-100%).  [required]")
	parser.add_option("-M","--maxFreq",
					  dest = "maxFreq",
					  default = "100",
					  metavar = "float",
					  help = "max Freq of snv to calculate (0-100%).  [required]")
	parser.add_option("-s","--step",
					  dest = "step",
					  default = "5",
					  metavar = "int",
					  help = "frame step to calculate.  [required]")

	parser.add_option("-T","--seqTech",
					  dest = "seqTech",
					  default = "",
					  metavar = "str",
					  help = "seq Tech (NGS CCS).  [required]")

	(options,args) = parser.parse_args()
	iSNV_SNP_table   = os.path.abspath(options.iSNV_SNP_table)
	out_P		  = os.path.abspath(options.out_P)
	RepeatRegionP = os.path.abspath(options.RepeatRegionP)
	minFreq		= float(options.minFreq)
	maxFreq		= float(options.maxFreq)
	step	    = int(options.step)
	seqTech     = options.seqTech

	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




