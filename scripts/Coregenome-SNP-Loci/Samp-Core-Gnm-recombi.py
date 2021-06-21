
#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys,os
from optparse import OptionParser


def SampCoreGnm(CoreGnmF,out_P):
	seqDic = {}
	for CoreGnmFl in open(CoreGnmF,'r').readlines():
		if ">" in CoreGnmFl and CoreGnmFl != "\n" :
			Samp = CoreGnmFl.split("\n")[0].split(">")[1]
			seqDic[Samp] = ''
			print(Samp)
		else:
			seqDic[Samp] += CoreGnmFl.split("\n")[0]	
	return seqDic

				
def allSamps_recombiRgs(ClonalFrameF):
	SampsRcombRgs = {}
	for ClonalFrameFl in open(ClonalFrameF).readlines():
		if "NODE_" not in ClonalFrameFl and "Node" not in ClonalFrameFl :
			Tags = ClonalFrameFl.split("\n")[0].split("\t")
			Samp = Tags[0]
			Beg = Tags[1]
			End = Tags[2]
			if Samp not in SampsRcombRgs:
				SampsRcombRgs[Samp] = []
			SampsRcombRgs[Samp].append(Beg + "_" + End)

	return SampsRcombRgs



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




def main():
	seqDic = SampCoreGnm(CoreGnmF,out_P)
	#personConsensusCoreGnm(CoreGnmF)
	#PersonCoreGnm_Diff(CoreGnmF)
	print(seqDic.keys())
	#SampsRcombRgs = allSamps_recombiRgs(ClonalFrameF)
	#print(SampsRcombRgs)


	personCont = {}
	for Sam in seqDic.keys():
		print(Sam)
		SEQ = seqDic[Sam]
		person = Sam.split("-")[0]
		if Sam not in  ["GCF_000224555.1_ASM22455v1","P24-C_sC-j","P8-E-t"]:
			if person not in personCont:
				personCont[person] = 1
		
			outF = out_P + "/"+ person + ".coreGnm-Seqs.fasta"
			if (personCont[person] == 1 and  os.path.exists(outF)):
				os.remove(outF)
			out_P_O=open(outF,'a')
			out_P_O.write(">" + Sam + "\n")
			out_P_O.write(SEQ + "\n")
			out_P_O.close()	
			personCont[person] +=1


'''
	personCont = {}
	persSampCont = {}
	for Sam in seqDic.keys():
		print(Sam)
		SEQ = seqDic[Sam]
		person = Sam.split("-")[0]
		if Sam not in  ["GCF_000224555.1_ASM22455v1","P24-C_sC-j","P8-E-t"]:
			if person not in personCont:
				personCont[person] = 0
			personCont[person] +=1
	
	count = {}
	for Sam in seqDic.keys():
		print(Sam)
		SEQ = seqDic[Sam]
		person = Sam.split("-")[0]
		if Sam not in  ["GCF_000224555.1_ASM22455v1","P24-C_sC-j","P8-E-t"]:
			ref_samp = PersonRespectSamp(person)
			if person not in count:
				count[person] = 0
			count[person] +=1

			
			outF = out_P + "/"+ person + ".coreGnm-Seqs.fasta"
			if (count[person] == 1 and  os.path.exists(outF)):
				os.remove(outF)
			print(count[person])
			print(personCont[person])
			if Sam != ref_samp:
				out_P_O=open(outF,'a')
				out_P_O.write(">" + Sam + "\n")
				out_P_O.write(SEQ + "\n")
				out_P_O.close()	
			else :
				for time in range(0,5-personCont[person]):
					print(time)
					out_P_O=open(outF,'a')
					out_P_O.write(">" + ref_samp + "__" + str(time) + "\n")
					out_P_O.write(seqDic[ref_samp] + "\n")
					out_P_O.close()	
					


	for Sam in seqDic.keys():
		print(Sam)
		SEQ = seqDic[Sam]
		outF = out_P + "/"+ Sam + "coreGnm-recombi-Seqs.fasta"
		if (os.path.exists(outF)):
			os.remove(outF)
		out_P_O=open(outF,'a')
		if Sam in SampsRcombRgs:			
			recombs = SampsRcombRgs[Sam]
			for recomb in recombs:
				Beg = int(recomb.split("_")[0])
				End = int(recomb.split("_")[1])
				recbSeq = SEQ[Beg-1:End]
				print(Sam)
				print(recomb)
				print(len(recbSeq))

				out_P_O.write(">" + Sam +  "--" +recomb + "\n")
				out_P_O.write(recbSeq + "\n")
		out_P_O.close()	

'''



	





if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-c","--CoreGnmF",
					  dest = "CoreGnmF",
					  default = "",
					  metavar = "file",
					  help = "core genome file from roary.  [required]")

	parser.add_option("-C","--ClonalFrameF",
					  dest = "ClonalFrameF",
					  default = "",
					  metavar = "file",
					  help = "ClonalFrame output files.  [required]")
	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output path.  [required]")



	(options,args) = parser.parse_args()
	CoreGnmF   = os.path.abspath(options.CoreGnmF)
	ClonalFrameF  = os.path.abspath(options.ClonalFrameF)
	out_P		  = os.path.abspath(options.out_P)



	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()





