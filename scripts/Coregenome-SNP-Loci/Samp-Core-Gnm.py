
#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys,os
from optparse import OptionParser
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


def SampCoreGnm(CoreGnmF,out_P):
	seqDic = {}
	for CoreGnmFl in open(CoreGnmF,'r').readlines():
		if ">" in CoreGnmFl and CoreGnmFl != "\n" :
			Samp = CoreGnmFl.split("\n")[0]
			seqDic[Samp] = ''
			print(Samp)
		else:
			seqDic[Samp] += CoreGnmFl.split("\n")[0]

	for SAMP in seqDic:
		outF = out_P + "/" + SAMP.split(">")[1] + "_CoreGnm.fasta"
		if ( os.path.exists(outF)):
			os.remove(outF)
		out_P_O=open(outF,'a')
		out_P_O.write(SAMP + "\n")
		out_P_O.write(seqDic[SAMP] + "\n")
		out_P_O.close()	

	return seqDic


def vcf_Header(SAMPID):
	HeaderL = []
	HeaderL.append("##fileformat=VCFv4.1")
	HeaderL.append("##source=iSNV")
	HeaderL.append('##INFO=<ID=ADP,Number=1,Type=Integer,Description="Average per-sample depth of bases with Phred score >= 20">')
	HeaderL.append('##INFO=<ID=WT,Number=1,Type=Integer,Description="Number of samples called reference (wild-type)">')
	HeaderL.append('##INFO=<ID=HET,Number=1,Type=Integer,Description="Number of samples called heterozygous-variant">')
	HeaderL.append('##INFO=<ID=HOM,Number=1,Type=Integer,Description="Number of samples called homozygous-variant">')
	HeaderL.append('##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of samples not called">')
	HeaderL.append(str('##FILTER=<ID=str10,Description="Less than 10% or more than 90% of variant supporting reads on one strand">'))
	HeaderL.append('##FILTER=<ID=indelError,Description="Likely artifact due to indel reads at this position">')
	HeaderL.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
	HeaderL.append('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
	HeaderL.append('##FORMAT=<ID=MDP,Number=1,Type=Integer,Description="Raw Read Depth from mpileup">')
	HeaderL.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= MapQThreds">')
	HeaderL.append('##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">')
	HeaderL.append('##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">')
	HeaderL.append('##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">')
	HeaderL.append('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	' + SAMPID)

	HEADER = "\n".join(HeaderL)
	return HEADER



def PersonCoreGnm_Diff(CoreGnmF):
	seqDic,personDic = {},{}
	
	for CoreGnmFl in open(CoreGnmF,'r').readlines():
		if ">" in CoreGnmFl and CoreGnmFl != "\n" :
			Samp = CoreGnmFl.split("\n")[0].split(">")[1]
			seqDic[Samp] = ''
			print(Samp)
			person = Samp.split("-")[0]
			if person not in personDic:
				personDic[person] = {}
			personDic[person][Samp] = ''
			seqLen = 0
		else:
			seqDic[Samp] += CoreGnmFl.split("\n")[0]
			seqLen += len(CoreGnmFl.split("\n")[0])
	print(seqLen)
		
	for person in personDic:
		if person != "GCF_000224555.1_ASM22455v1":
			ref_samp = PersonRespectSamp(person)
			ref_sampGnm = seqDic[ref_samp]
			print(person)
			print(ref_samp)
			for Samp in personDic[person]:
				outL = []
				#if Samp != ref_samp:
				HEADER = vcf_Header(Samp)
				for Posinx in range(0,seqLen):
					refBase = ref_sampGnm[Posinx]
					SampBs = seqDic[Samp][Posinx]
					if SampBs != refBase and refBase != "-" and SampBs!= "-":
						outlineL = []
						outlineL.append(ref_samp)
						outlineL.append(str(Posinx+1))
						outlineL.append(".")
						outlineL.append(refBase)  ##REF
						outlineL.append(SampBs)  ##ALT
						outlineL.append(".")
						outlineL.append("PASS")
						outlineL.append(";".join(["ADP=0","WT=0","HET=0","HOM=1","NC=0"]))
						outlineL.append("GT:GQ:MDP:DP:RD:AD:FREQ")
						outlineL.append(":".join(["1/1","255","0","0","0","0","0"]))							
						outline = "\t".join(outlineL)
						outL.append(outline)
				outF = out_P + "/" + Samp + "_CoreGnm_Diffposi_map2Respect.vcf"
				if ( os.path.exists(outF)):
					os.remove(outF)
				out_P_O=open(outF,'a')
				out_P_O.write(HEADER + "\n")
				out_P_O.write("\n".join(outL) + "\n")
				out_P_O.close()	



'''
def personConsensusCoreGnm(CoreGnmF):
	seqDic,personDic = {},{}
	
	for CoreGnmFl in open(CoreGnmF,'r').readlines():
		if ">" in CoreGnmFl and CoreGnmFl != "\n" :
			Samp = CoreGnmFl.split("\n")[0].split(">")[1]
			seqDic[Samp] = ''
			print(Samp)
			person = Samp.split("-")[0]
			if person not in personDic:
				personDic[person] = {}
			personDic[person][Samp] = ''
			seqLen = 0
		else:
			seqDic[Samp] += CoreGnmFl.split("\n")[0]
			seqLen += len(CoreGnmFl.split("\n")[0])
	print(seqLen)
	consensusSeqDic = {}		
	for person in personDic:
		print(person)
		consensusSeqDic[person] = ''
		for Posinx in range(0,seqLen):
			BaseDic = {}
			for Samp in personDic[person]:
				SampBs = seqDic[Samp][Posinx]
				if SampBs not in BaseDic:
					BaseDic[SampBs] = 0
				BaseDic[SampBs] += 1
			Maxbase = max(BaseDic,key=BaseDic.get)
			#if len(BaseDic) > 1:
				#print(BaseDic)
				#print(Maxbase)
			consensusSeqDic[person] += Maxbase
	for person in personDic:
		outF = out_P + "/" + person + "_CoreGnm_consensusSeq.fasta"
		if ( os.path.exists(outF)):
			os.remove(outF)
		out_P_O=open(outF,'a')
		out_P_O.write(">" + person + "\n")
		out_P_O.write(consensusSeqDic[person] + "\n")
		out_P_O.close()	

'''
				





def main():
	SampCoreGnm(CoreGnmF,out_P)
	#personConsensusCoreGnm(CoreGnmF)
	PersonCoreGnm_Diff(CoreGnmF)




if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-c","--CoreGnmF",
					  dest = "CoreGnmF",
					  default = "",
					  metavar = "file",
					  help = "core genome file from roary.  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output path.  [required]")



	(options,args) = parser.parse_args()
	CoreGnmF   = os.path.abspath(options.CoreGnmF)
	out_P		  = os.path.abspath(options.out_P)



	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()





