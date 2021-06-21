
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
		if lie > 14 and  Header != "P11-E-j.ONT" and Header != "P8-E-t" and Header != "P24-C_sC-j":			 
			#Person = Header.split("-")[0]
			#if Person not in PersonDic:
				#PersonDic[Person] = {}
			PersonDic[Header] = ''
	print(PersonDic)

	Genegrop = {}
	for pangenomel in pangenomels:		
		if "Gene" not in pangenomel and pangenomel != "\n" :
			lL = pangenomel.split("\n")[0].split("\",\"")
			genegrpID = lL[0].split("\"")[1]	
			GeneLSt = {}
			for INdex in range(14,len(lL)):
				sampID = lieD[INdex]
				if sampID != "P11-E-j.ONT" and sampID != "P8-E-t" and sampID != "P24-C_sC-j":			 
					Gid = lL[INdex]
					if INdex == len(lL)-1:
						Gid = lL[INdex].split("\"")[0]
					if Gid != '':
						if Gid not in GeneLSt :
							GeneLSt[sampID] = Gid
			if len(GeneLSt) > 1:
				Genegrop[genegrpID] = GeneLSt
	return PersonDic,Genegrop


def readNucF(person,samp,prokkaFileP,personNucD):
	personNucD[samp] = {}
	NucFile = prokkaFileP + "/" + samp + "_prokka/"+ samp + ".ffn"
	print(NucFile)
	for l in open(NucFile).readlines():
		if l != "\n":
			if ">" in l:
				H = l.split("\n")[0].split(" ")[0].split(">")[1]
				personNucD[samp][H] = ''
			else:
				personNucD[samp][H] += l.split("\n")[0]
	return personNucD


	
def readAAF(person,samp,AAFileP,personAAD):
	personAAD[samp] = {}
	AAFile = AAFileP + "/" + samp + "_prokka/" + samp + ".faa"
	print(AAFile)
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

	

	personNucD,personAAD = {},{}
	for samp in PersonDic:
		print(samp)
		sampProkkaP = prokkaFileP
		readNucF(personID,samp,sampProkkaP,personNucD)
		readAAF(personID,samp,sampProkkaP,personAAD)

	#print(personNucD)
	for GenegropID in Genegrop:
		GrpGens = []
		GrpAAs = []
		personGeneLst = list(Genegrop[GenegropID].keys())
		firstSamp = personGeneLst[0]
		#print(firstSamp)
		firstGne = Genegrop[GenegropID][firstSamp]
		#print(firstGne)
		firstGneSeq = personNucD[firstSamp][firstGne]
		diffCount = 0
		GrpGens.append(">" + GenegropID + "__" + firstSamp + "__" + firstGne)
		GrpGens.append(firstGneSeq)
		
		firstAASeq = personAAD[firstSamp][firstGne]	
		GrpAAs.append(">" + GenegropID + "__" + firstSamp + "__" + firstGne)
		GrpAAs.append(personAAD[firstSamp][firstGne])

		print(personAAD.keys())
		for OtherSamp in personGeneLst[1:]:
			#Nuc
			OtherSampGne = Genegrop[GenegropID][OtherSamp]
			OtherSampGneSeq = personNucD[OtherSamp][OtherSampGne]	
			GrpGens.append(">" + GenegropID + "__" + OtherSamp + "__" + OtherSampGne)
			GrpGens.append(OtherSampGneSeq)
			if OtherSampGneSeq != firstGneSeq:
				diffCount +=1
			
			##AA
			OtherSampAASeq = personAAD[OtherSamp][OtherSampGne]	
			GrpAAs.append(">" + GenegropID + "__" + OtherSamp + "__" + OtherSampGne)
			GrpAAs.append(OtherSampAASeq)


		if len(personGeneLst) >= 2 and diffCount != 0 :  #and "group_13" in GenegropID
			NucOut = "\n".join(GrpGens)
			print(NucOut)

			outNucF = out_P+ "/" + personID  + "__" + GenegropID + ".nuc.fna"
			if ( os.path.exists(outNucF)):
				os.remove(outNucF)
			out_P_O=open(outNucF,'a')
			out_P_O.write(NucOut + "\n")
			out_P_O.close()	


			AAOut = "\n".join(GrpAAs)
			print(NucOut)

			AAoutF = out_P+ "/" + personID  + "__" + GenegropID + ".aa.fna"
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





