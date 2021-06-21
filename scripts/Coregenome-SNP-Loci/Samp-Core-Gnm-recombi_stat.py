
#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys,os
from optparse import OptionParser


def SampCoreGnm(CoreGnmF):
	seqDic = {}
	count = 0
	for CoreGnmFl in open(CoreGnmF,'r').readlines():
		if ">" in CoreGnmFl and CoreGnmFl != "\n" :
			Samp = CoreGnmFl.split("\n")[0].split(">")[1]
			seqDic[Samp] = ''
			print(Samp)
			count += 1
			if count == 2 :
				break
		else:
			seqDic[Samp] += CoreGnmFl.split("\n")[0]
			Len = len(seqDic[Samp])	

	return Len

				
def allSamps_recombiRgs(ClonalFrameP,out_P,Corelen):
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
	Pls = []
	Pls.append("\t".join(["person","Corelen","RecombPosNum"]))
	for person in personRecombi:
		Pl = "\t".join([person,str(Corelen),str(len(personRecombi[person]))])
		Pls.append(Pl)
	PO = "\n".join(Pls)
	PoutF = out_P + "/person.CoreGnm.Stat.txt"
	if (os.path.exists(PoutF)):
		os.remove(PoutF)
	out_P_O=open(PoutF,'a')
	out_P_O.write(PO + "\n")
	out_P_O.close()	



	Sls = []
	Sls.append("\t".join(["person","Corelen","RecombPosNum"]))
	for Samp in SampRecombi:
		Sl = "\t".join([Samp,str(Corelen),str(len(SampRecombi[Samp]))])
		Sls.append(Sl)
	SO = "\n".join(Sls)

	SoutF = out_P + "/Samps.CoreGnm.Stat.txt"
	if (os.path.exists(SoutF)):
		os.remove(SoutF)
	out_P_O=open(SoutF,'a')
	out_P_O.write(SO + "\n")
	out_P_O.close()	



def main():
	Corelen = SampCoreGnm(CoreGnmF)

	allSamps_recombiRgs(ClonalFrameP,out_P,Corelen)
	





if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-c","--CoreGnmF",
					  dest = "CoreGnmF",
					  default = "",
					  metavar = "file",
					  help = "core genome file from roary.  [required]")

	parser.add_option("-C","--ClonalFrameP",
					  dest = "ClonalFrameP",
					  default = "",
					  metavar = "path",
					  help = "ClonalFrame output files.  [required]")
	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output path.  [required]")



	(options,args) = parser.parse_args()
	CoreGnmF   = os.path.abspath(options.CoreGnmF)
	ClonalFrameP  = os.path.abspath(options.ClonalFrameP)
	out_P		  = os.path.abspath(options.out_P)



	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()





