#!/usr/bin/python
# -*- coding: utf-8 -*-

# Complementing DNA 
import sys,os
import linecache
from optparse import OptionParser



def readVCF(sample,SV_vcfF,SV_Dic,allSV_dic):
	SV_Dic[sample] = {}
	allSV_dic[sample] = {}
	SV_vcfFls=open(SV_vcfF).readlines()
	for SV_vcfFl in SV_vcfFls:
		if "#" not in SV_vcfFl and SV_vcfFl != "\n":
			lTags = SV_vcfFl.split("\t")
			Posi = lTags[1]
			#SVType=lTags[4]
			Infos = lTags[7].split(";")
			
			for Info in Infos:
				if "SVLEN" in Info:
					SV_len = Info.split("=")[-1]
				if "SVTYPE" in Info:
					SVType = Info.split("=")[-1]

			if abs(int(SV_len)) < maxLen :
				if "DEL" in SVType:
					AlleDP=int(lTags[9].split("\n")[0].split(":")[-1])
					if AlleDP >= minADP:
						if Posi not in SV_Dic[sample]:
							SV_Dic[sample][Posi] = {}
						SV_Dic[sample][Posi]["DEL"] = int(SV_len)

						allSV_dic[sample][Posi] = SV_vcfFl

				elif "INS" in SVType:
					AlleDP=int(lTags[9].split("\n")[0].split(":")[-1])
					if AlleDP >= minADP:
						if Posi not in SV_Dic[sample]:
							SV_Dic[sample][Posi] = {}
						SV_Dic[sample][Posi]["INS"] = int(SV_len)

						allSV_dic[sample][Posi] = SV_vcfFl

	return SV_Dic,allSV_dic	




def main():
	personDic={}
	vcfFiles = os.listdir(vcfP)
	countLs = []
	#print(vcfFiles)
	for file in vcfFiles:
		sampID = file.split("/")[-1].split(".")[0]
		person = sampID.split("-")[0]
		if person not in personDic:
			personDic[person] = []
		if "P24-C_sC-j" not in  sampID :
			personDic[person].append(sampID)

	for Person in personDic:
		SV_Dic,allSV_dic,person_SV_Cites = {},{},{}
		SV_Cites_count = {}
		head = ["Cites","sharedNum"]

		for sample in personDic[Person]:
			head.append(sample)
			for file in vcfFiles:
				if sample in file :
					SV_vcfF = vcfP +  "/" +file
			SV_Dic,allSV_dic = readVCF(sample,SV_vcfF,SV_Dic,allSV_dic)
			
			#print(SV_Dic[sample])
			for Posi in SV_Dic[sample]:
				if Posi not in person_SV_Cites:
					person_SV_Cites[Posi] = ''
				if Posi not in SV_Cites_count:
					SV_Cites_count[Posi] = 0
				SV_Cites_count[Posi] += 1
		outls,SV_lines = [],[]	


		#print()

		DelCount,INSCount = {},{}

		for Cite in person_SV_Cites:
			outl = [Cite,str(SV_Cites_count[Cite])]
			for sample in personDic[Person]:
				#if SV_Cites_count[Cite] != len(personDic[Person]) :			
				if sample not in DelCount:
					DelCount[sample] = 0
				if sample not in INSCount:
					INSCount[sample] = 0
				#if SV_Cites_count[Cite] != len(personDic[Person]) :
				if Cite in SV_Dic[sample]:
					if "DEL" in SV_Dic[sample][Cite]:
						outl.append("DEL." + str(SV_Dic[sample][Cite]["DEL"]) + "(" + str(allSV_dic[sample][Cite].split("\t")[9].split("\n")[0].split(":")[-1]) + ")")
						SVtype="DEL"
						SVLen = int(SV_Dic[sample][Cite]["DEL"])
						#if sample not in DelCount:
							#DelCount[sample] = 0
						DelCount[sample] += 1
					elif "INS" in SV_Dic[sample][Cite]:
						outl.append("INS." + str(SV_Dic[sample][Cite]["INS"]) + "(" + str(allSV_dic[sample][Cite].split("\t")[9].split("\n")[0].split(":")[-1]) + ")" )
						SVtype="INS"
						SVLen = int(SV_Dic[sample][Cite]["INS"])						
						#if sample not in INSCount:
							#INSCount[sample] = 0
						INSCount[sample] += 1
						print(INSCount)
					SV_lines.append(sample + "\t" + SVtype + "\t" + str(SVLen) + "\t" + allSV_dic[sample][Cite])
				else:
					outl.append('-')
			if len(outl) != 1:
				outls.append("\t".join(outl))

		#print(DelCount)
		for sample in personDic[Person]:
			countLs.append("\t".join([sample.split("_2_")[0],sample.split("_2_")[1],str(DelCount[sample]),str(INSCount[sample] )]))

		#print("\t".join(head))
		#print("\n".join(outls))
		out_F = out_P + "/" + Person + ".CCS-SV.sniffles.minADP-" + str(minADP) + ".DEL-INS.txt"

		if (os.path.exists(out_F)):
			os.remove(out_F)

		out_F_O=open(out_F,'a')
		out_F_O.write("\t".join(head) + "\n")	
		out_F_O.write("\n".join(outls)+ "\n")	
		out_F_O.close()	


		out_filtedSV_F = out_P + "/" + Person + ".CCS-SV.sniffles.minADP-" + str(minADP) + ".DEL-INS.filtedSV.txt"
		if (os.path.exists(out_filtedSV_F)):
			os.remove(out_filtedSV_F)

		out_F_O=open(out_filtedSV_F,'a')
		out_F_O.write("".join(SV_lines)	)
		out_F_O.close()	



		out_count_F = out_P + "/allSamps.DEL-INS.count.txt"
		if (os.path.exists(out_count_F)):
			os.remove(out_count_F)
		out_F_O=open(out_count_F,'a')
		out_F_O.write("\n".join(countLs)	)
		out_F_O.close()	





if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-v","--vcfP",
					  dest = "vcfP",
					  default = "",
					  metavar = "path",
					  help = "path of vcf files.  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output stat file .  [required]")

	parser.add_option("-m","--minADP",
					  dest = "minADP",
					  default = "50",
					  metavar = "int",
					  help = "min depth support the alle.  [required]")

	parser.add_option("-M","--maxLen",
					  dest = "maxLen",
					  default = "10000",
					  metavar = "int",
					  help = "max length of the alle  [required]")

	(options,args) = parser.parse_args()
	vcfP   = os.path.abspath(options.vcfP)
	out_P		  = os.path.abspath(options.out_P)
	minADP		= int(options.minADP)
	maxLen		= int(options.maxLen)
	
	
	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()










