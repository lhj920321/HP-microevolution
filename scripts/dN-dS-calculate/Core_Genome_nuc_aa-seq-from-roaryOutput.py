
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 


def Read_pan_genes(pan_genes_csv):
	pangenomels = open(pan_genes_csv).readlines()
	HeaderL = pangenomels[0].split("\n")[0].split("\",\"")
	lie = 0
	PersonDic = {}
	lieD = {}
	for Header in HeaderL:
		if Header in HeaderL[-1]:
			Header= Header.split("\"")[0]
		lieD[lie] = Header
		lie += 1

		#if lie > 14 and  Header != "P11-E-j.ONT" and Header != "P8-E-t" and Header != "P24-C_sC-j":			 
		if lie > 14 and Header != "P8-E-t" and Header != "P24-C_sC-j":			 
			Person = Header.split("-")[0]
			if Person not in PersonDic:
				PersonDic[Person] = {}
			PersonDic[Person][Header] = ''



	Genegrop = {}
	for pangenomel in pangenomels:		
		if "Gene" not in pangenomel and pangenomel != "\n" :
			lL = pangenomel.split("\n")[0].split("\",\"")
			genegrpID = lL[0].split("\"")[1]	
			GeneLSt = {}
			for INdex in range(14,len(lL)):
				sampID = lieD[INdex]
				if sampID != "P8-E-t" and sampID != "P24-C_sC-j":	  #sampID != "P11-E-j.ONT" and 		 
					Gid = lL[INdex]				
					if INdex == len(lL)-1:
						Gid = lL[INdex].split("\"")[0]			
					if Gid != '':
						if sampID not in GeneLSt :
							GeneLSt[sampID] = Gid
		
			if len(GeneLSt) >= 70:
				Genegrop[genegrpID] = GeneLSt

	return PersonDic,Genegrop





def readNucF(samp,prokkaFileP,personNucD):
	personNucD[samp] = {}
	NucFile = prokkaFileP + "/" + samp + "_prokka/"+ samp + ".ffn"
	#print(NucFile)
	for l in open(NucFile).readlines():
		if l != "\n":
			if ">" in l:
				H = l.split("\n")[0].split(" ")[0].split(">")[1]
				personNucD[samp][H] = ''
			else:
				personNucD[samp][H] += l.split("\n")[0]
	return personNucD


	
def readAAF(samp,AAFileP,personAAD):
	personAAD[samp] = {}
	AAFile = AAFileP + "/" + samp + "_prokka/" + samp + ".faa"
	for l in open(AAFile).readlines():
		if l != "\n":
			if ">" in l:
				H=l.split("\n")[0].split(" ")[0].split(">")[1]
				personAAD[samp][H] = ''
			else:
				personAAD[samp][H] += l.split("\n")[0]
	return personAAD



					





def main():
	pan_genes_csv = roaryP+ "/gene_presence_absence.csv"
	PersonDic,Genegrop = Read_pan_genes(pan_genes_csv)

	
	#print(Genegrop)
	
	for personIDD in PersonDic:
		print(personIDD)
		if personIDD == personID:
			personNucD,personAAD = {},{}
			if  "GCF_000224555.1_ASM22455v1" not in personIDD :
				for samp in PersonDic[personIDD]:
					sampProkkaP = prokkaFileP
					personNucD=readNucF(samp,sampProkkaP,personNucD)
					personAAD = readAAF(samp,sampProkkaP,personAAD)

				for GenegropID in Genegrop:
					GrpGens = []
					GrpAAs = []
					for SampGene in Genegrop[GenegropID]:
						if  "GCF_000224555.1_ASM22455v1" not in SampGene :
							#Nuc
							SampGne = Genegrop[GenegropID][SampGene]
							if SampGene in personNucD and   SampGne in personNucD[SampGene]:
								SampGneSeq = personNucD[SampGene][SampGne]	
								GrpGens.append(">" + GenegropID + "__" + SampGene + "__" + SampGne)
								GrpGens.append(SampGneSeq)
								
							##AA
							if SampGene in personNucD and  SampGne in personAAD[SampGene]:
								SampAASeq = personAAD[SampGene][SampGne]	
								GrpAAs.append(">" + GenegropID + "__" + SampGene + "__" + SampGne)
								GrpAAs.append(SampAASeq)


					NucOut = "\n".join(GrpGens)
					#OP = "/".join(out_P.split("/")[:-1]) + "/CoreGnm_" + personIDD 
					#if ( not os.path.exists(OP)):
						#os.mkdir(OP)

					outNucF = out_P + "/" + personID  + "__" + GenegropID + ".nuc.fna"
					#print(outNucF)
					if ( os.path.exists(outNucF)):
						os.remove(outNucF)
					out_P_O=open(outNucF,'a')
					out_P_O.write(NucOut + "\n")
					out_P_O.close()	
					AAOut = "\n".join(GrpAAs)
					#print(NucOut)
					AAoutF = out_P + "/" + personID  + "__" + GenegropID + ".aa.fna"
					if ( os.path.exists(AAoutF)):
						os.remove(AAoutF)
					out_P_O=open(AAoutF,'a')
					out_P_O.write(AAOut + "\n")
					out_P_O.close()	
					




if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--personID",
					  dest = "personID",
					  default = "",
					  metavar = "str",
					  help = "person ID .  [required]")


	parser.add_option("-r","--roaryP",
					  dest = "roaryP",
					  default = "",
					  metavar = "path",
					  help = "Path of roary output files .  [required]")

	parser.add_option("-p","--prokkaFileP",
					  dest = "prokkaFileP",
					  default = "",
					  metavar = "path",
					  help = "Path of sample's prokka output files  .  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output stat file Path of snv Freq distribution.  [required]")



	(options,args) = parser.parse_args()
	personID   = options.personID
	roaryP   = os.path.abspath(options.roaryP)
	prokkaFileP   = os.path.abspath(options.prokkaFileP)
	out_P		  = os.path.abspath(options.out_P)



	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)



	main()





