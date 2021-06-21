
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 


def PersonAnnoPosi(GffP):
	PersonDiffPoAnno = {}
	for file in os.listdir(GffP):
		#print(file)
		Samp = file.split("_2_")[0]
		person = Samp.split("-")[0]
		if Samp.split("_2_")[0] not in ["GCF_000224555.1_ASM22455v1","P24-C_sC-j","P8-E-t"]:
			if person not in PersonDiffPoAnno:
				PersonDiffPoAnno[person] = {}
			for Fl in open(GffP+"/"+file).readlines() :
				if "#" not in Fl  and Fl != "\n":
					Posi = Fl.split("\t")[1]
					Anno = Fl.split("\t")[7].split(";")[5].split("|")[1]
					if Posi not in PersonDiffPoAnno[person]:
						PersonDiffPoAnno[person][Posi] = Anno
					else:
						if Anno != PersonDiffPoAnno[person][Posi]:
							print("nnonnononono")
							print(Posi)
							print(Anno)
							print(PersonDiffPoAnno[person][Posi])
							print()
	return PersonDiffPoAnno

'''
def read_Recomb(RecombP):
	RecombPosis = {}
	for file in os.listdir(RecombP):
		Samp = file.split(".")[0]
		person = Samp.split("-")[0]
		RecombPosis[Samp] = {}
		for Fl in open(GffP+"/"+file).readlines() :
			if "Node" not in Fl and ("NODE" not in Fl) and Fl != "\n":
				Tags = ClonalFrameFl.split("\n")[0].split("\t")
				Beg = int(Tags[1])
				End = int(Tags[2])
				for CT in range(Beg,End + 1):
					RecombPosis[Samp][CT] = ''


	persRecombiPosi = {}		
	
		SampRecombPosis = RecombPosis[Samp]   ##RecombPosis
		for cite in SampRecombPosis:
			if cite not in persRecombiPosi:
				persRecombiPosi[cite] = ''

'''


		


def Gnms_stat(CoreGnmF,PersonDiffPoAnno,SampRecombi,personRecombi):
	PerSeqsdic = {}
	for CoreGnmFl in open(CoreGnmF,'r').readlines():
		if ">" in CoreGnmFl and CoreGnmFl != "\n" :
			Samp = CoreGnmFl.split("\n")[0].split(">")[1].split(".pilon.3")[0]
			print(Samp)
			if Samp not in ["GCF_000224555.1_ASM22455v1","P24-C_sC-j","P8-E-t"]:
				Person = Samp.split("-")[0]
				if Person not in PerSeqsdic:
					PerSeqsdic[Person] = {}
				if Samp not in PerSeqsdic[Person]:
					PerSeqsdic[Person][Samp] = ''

		else:
			if Samp not in ["GCF_000224555.1_ASM22455v1","P24-C_sC-j","P8-E-t"]:
				PerSeqsdic[Person][Samp] += CoreGnmFl.split("\n")[0]
				SeqLen = len(PerSeqsdic[Person][Samp])
	print(SeqLen)

	print("core Seq Len : " + str(SeqLen))
	personNAPosiNum,personDiffPosiNum,personDiffAnnoDic,RecombiDiffNum = {},{},{},{}
	for personID in PerSeqsdic:
		personDiffPosiNum[personID],personNAPosiNum[personID]  = 0,0
		personDiffAnnoDic[personID] = {}
		print(personID)
		RecombiDiffNum[personID] = 0
		 
		for CorePos in range(0,SeqLen):
			
			Pers_Bases = {}
			for Samp in PerSeqsdic[personID]:
				sampBase = PerSeqsdic[personID][Samp][CorePos]
				if sampBase not in Pers_Bases:
					Pers_Bases[sampBase] = 0
				Pers_Bases[sampBase] += 1
			if "-" not in Pers_Bases and len(Pers_Bases) != 1:
				personDiffPosiNum[personID] += 1
				#print(Pers_Bases)
				Anno = PersonDiffPoAnno[personID][str(CorePos+1)]
				if Anno not in personDiffAnnoDic[personID]:
					personDiffAnnoDic[personID][Anno] = 0
				personDiffAnnoDic[personID][Anno] += 1
				if CorePos+1 in personRecombi[personID]:
					RecombiDiffNum[personID] += 1					
			if "-" in Pers_Bases and len(Pers_Bases) != len(PerSeqsdic[personID]):
				personNAPosiNum[personID] += 1
	
	for pers in PerSeqsdic:
		print(pers+ "\t"+ str(personNAPosiNum[pers]))
	print()
	print()
	print(personDiffAnnoDic[personID].keys())
	print()
	outs,outAnno = [],[]
	for pers in PerSeqsdic:
		if int(pers.split("P")[1]) >= 1 and int(pers.split("P")[1]) <= 14:
			group = "E"
		else:
			group = "C"
	
		SynNum,MissNum,StopGainNum,StopLostNum,otherAnno  = 0,0,0,0,0
		for anno in personDiffAnnoDic[pers]:
			if "synonymous_variant" in anno :
				SynNum += personDiffAnnoDic[pers][anno] 
			elif "missense_variant" in anno :
				MissNum += personDiffAnnoDic[pers][anno] 
			elif "stop_gained" in anno :
				StopGainNum += personDiffAnnoDic[pers][anno] 
			elif "stop_lost" in anno :
				StopLostNum += personDiffAnnoDic[pers][anno] 
			else :
				otherAnno +=  personDiffAnnoDic[pers][anno] 

		outl = "\t".join([pers,group,str(personDiffPosiNum[pers]),str(StopGainNum),str(StopLostNum),str(otherAnno),str(SynNum),str(MissNum),str(RecombiDiffNum[pers])])
		outs.append(outl)
		outAnnol = "\t".join([pers,str(StopGainNum),str(StopLostNum),str(otherAnno),str(SynNum),str(MissNum)])
		outAnno.append(outAnnol)
				
	return 	outs,outAnno

				
def allSamps_recombiRgs(ClonalFrameP):
	SampRecombi,personRecombi = {},{}
	for file in os.listdir(ClonalFrameP):
		if "importation_status.txt" in file:	
			print(file)
			person = file.split(".")[1]
			ClonalFrameF = ClonalFrameP + "/" + file
				
			if person not in personRecombi:
				personRecombi[person] = {}
			for ClonalFrameFl in open(ClonalFrameF).readlines():
				if "NODE_" not in ClonalFrameFl and "Node" not in ClonalFrameFl :
					Tags = ClonalFrameFl.split("\n")[0].split("\t")
					Samp = Tags[0]
					Beg = int(Tags[1])
					End = int(Tags[2])
					if Samp not in SampRecombi:
						SampRecombi[Samp] = {}
					for Pos in range(Beg,End):
						if Pos not in SampRecombi[Samp]:
							SampRecombi[Samp][Pos] = ''
						if Pos not in personRecombi[person]:
							personRecombi[person][Pos] = ''
	return SampRecombi,personRecombi


	

def main():
	PersonDiffPoAnno = PersonAnnoPosi(GffP)

	SampRecombi,personRecombi = allSamps_recombiRgs(ClonalFrameP)

	outs,outAnno = Gnms_stat(CoreGnmF,PersonDiffPoAnno,SampRecombi,personRecombi)
	OUts =[]
	OUts.append("\t".join(["person","group","SNPLoci","StopGainNum","StopLostNum","otherAnno","SynNum","MissNum","DiffposInRecomb"]))
	OUts += outs
	out_F_O=open(out_F,'a')
	out_F_O.write("\n".join(OUts) + "\n")
	out_F_O.close()	

	
	OutAnnos =[]
	OutAnnos.append("\t".join(["person","StopGainNum","StopLostNum","otherAnno","SynNum","MissNum"]))
	OutAnnos += outAnno
	O = out_F.split(".txt")[0] + ".Anno.txt"
	if ( os.path.exists(O)):
		os.remove(O)
	out_F_O=open(O,'a')
	out_F_O.write("\n".join(OutAnnos) + "\n")
	out_F_O.close()	



	H = "\t".join(["person","Anno","AnnoNum"])
	CO = out_F.split(".txt")[0]+ ".Cgroup.Anno.txt"
	if ( os.path.exists(CO)):
		os.remove(CO)
	C_O=open(CO,'a')
	C_O.write(H + "\n")

	EO = out_F.split(".txt")[0]+ ".Egroup.Anno.txt"
	if ( os.path.exists(EO)):
		os.remove(EO)
	E_O=open(EO,'a')
	E_O.write(H + "\n")

	for Ol in outAnno:
		if int(Ol.split("\t")[0].split("P")[1]) < 15:
			C_O.write( "\t".join([Ol.split("\t")[0],"StopGain",Ol.split("\t")[1]])  + "\n")
			C_O.write( "\t".join([Ol.split("\t")[0],"StopLostNum",Ol.split("\t")[2]])  + "\n")
			C_O.write( "\t".join([Ol.split("\t")[0],"otherAnno",Ol.split("\t")[3]])  + "\n")
			C_O.write( "\t".join([Ol.split("\t")[0],"SynNum",Ol.split("\t")[4]])  + "\n")
			C_O.write( "\t".join([Ol.split("\t")[0],"MissNum",Ol.split("\t")[5]])  + "\n")
		else:
			E_O.write( "\t".join([Ol.split("\t")[0],"StopGain",Ol.split("\t")[1]])  + "\n")
			E_O.write( "\t".join([Ol.split("\t")[0],"StopLostNum",Ol.split("\t")[2]])  + "\n")
			E_O.write( "\t".join([Ol.split("\t")[0],"otherAnno",Ol.split("\t")[3]])  + "\n")
			E_O.write( "\t".join([Ol.split("\t")[0],"SynNum",Ol.split("\t")[4]])  + "\n")
			E_O.write( "\t".join([Ol.split("\t")[0],"MissNum",Ol.split("\t")[5]])  + "\n")

	C_O.close()	
	E_O.close()	


if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-c","--CoreGnmF",
					  dest = "CoreGnmF",
					  default = "",
					  metavar = "file",
					  help = "core genome file from roary.  [required]")

	parser.add_option("-g","--GffP",
					  dest = "GffP",
					  default = "",
					  metavar = "path",
					  help = "Gff Path.  [required]")
	
	parser.add_option("-o","--out_F",
					  dest = "out_F",
					  default = "",
					  metavar = "file",
					  help = "output stat file Path of snv Freq distribution.  [required]")


	parser.add_option("-C","--ClonalFrameP",
					  dest = "ClonalFrameP",
					  default = "",
					  metavar = "path",
					  help = "ClonalFrame output files.  [required]")

	(options,args) = parser.parse_args()
	GffP  = os.path.abspath(options.GffP)
	CoreGnmF  = os.path.abspath(options.CoreGnmF)
	ClonalFrameP  = os.path.abspath(options.ClonalFrameP)
	out_F	  = os.path.abspath(options.out_F)



	if ( os.path.exists(out_F)):
		os.remove(out_F)


	main()





