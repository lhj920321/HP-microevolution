#!/usr/bin/python
# -*- coding: utf-8 -*-

###############
#modules
###############
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


def sample_ID(iSNVTable):
	##sample col  index  in table "iSNV_with_SNP.all.txt" 
	iSNVTablels = open(iSNVTable,'r').readlines()
	lieCount = 2
	for iSNVTablel in iSNVTablels:
		if "pos	" in iSNVTablel:
			sampIDs = iSNVTablel.split("\n")[0].split("\t")[2:-1]
			for sampID in sampIDs:
				lieCount += 1			
	return sampIDs  



def iSNVsampPos(iSNVTable,SAMPID,FREQ_threshold):
	iSNVTablels = open(iSNVTable,'r').readlines()
	##find sample col  index
	for iSNVTablel in iSNVTablels:
		if "pos	" in iSNVTablel:
			sampIDs = iSNVTablel.split("\n")[0].split("\t")
			sampIDcount = 0
			for sampID in sampIDs:
				if sampID == SAMPID:
					sampIDindex = sampIDcount
					break
				sampIDcount += 1
			break	

	###sample  posi  that isn't  NO or NA
	sampPosL = []
 	for iSNVTablel in iSNVTablels:
		if "#" not in iSNVTablel and iSNVTablel != "\n":
			sampsFreq = iSNVTablel.split("\n")[0].split("\t")
			POS = sampsFreq[0]
			sampFREQ = sampsFreq[sampIDindex]
			if sampFREQ != "NA" and sampFREQ != "NO":
				if float(sampFREQ) <= FREQ_threshold:
					sampPosL.append(POS)
	return sampPosL




def samp_posi_info(infoF,sampPosL,SAMPID,baseCovThreshold):
	infoFls = open(infoF,'r').readlines()
	sampSiteInfoD = {}
	notHegeL = []
	for infoFl in infoFls:
		if "#" not in infoFl and infoFl != "\n":
			infos = infoFl.split("\t")
			ID = infos[0]
			Site = infos[1]	

			if ID == SAMPID  and Site in sampPosL:	
				REF = infos[4].split(":")[0]
				nt_pattern = infos[18].split(";")
				Anum=int(nt_pattern[0].split(":")[1])
				Gnum=int(nt_pattern[1].split(":")[1])
				Cnum=int(nt_pattern[2].split(":")[1])
				Tnum=int(nt_pattern[3].split(":")[1])
				totnum=int(nt_pattern[4].split(":")[1])
				charsCount = 0
				NumD = {"A":Anum,"G":Gnum,"C":Cnum,"T":Tnum}
				FreqD = {}
				Afreq = str(round(float(Anum)*100/totnum,3))
				FreqD["A"] = Afreq
				Gfreq = str(round(float(Gnum)*100/totnum,3))
				FreqD["G"] = Gfreq
				Cfreq = str(round(float(Cnum)*100/totnum,3))
				FreqD["C"] = Cfreq
				Tfreq = str(round(float(Tnum)*100/totnum,3))
				FreqD["T"] = Tfreq
					
				if Anum > baseCovThreshold:	
					if REF != "A":
						charsCount += 1
				if Gnum > baseCovThreshold:
					if REF != "G":
						charsCount += 1
				if Cnum > baseCovThreshold:
					if REF != "C":
						charsCount += 1
				if Tnum > baseCovThreshold:
					if REF != "T":
						charsCount += 1


				if charsCount <= 3:   #max num of mutation type for one position 
				##if charsCount >= 1:   
					MDP = infos[5]   #cov_from_mpileup
					DP = infos[7]   #cov_noIndel_QC
					minor_cov = infos[9]
					second_cov = infos[11]	
					minor_freq = infos[10]
					second_freq = infos[12]						
					chars_max = infos[4].split(":")[1].split("-")[0]
					chars_second = infos[4].split(":")[1].split("-")[1]
					FREQ = infos[2]
					if chars_max == REF:
						ALT = chars_second						
					elif chars_second == REF:
						ALT = chars_max
					else:
						ALT = chars_max
					RD = str(NumD[REF])
					AD = str(NumD[ALT])
					#print RD + "\t"  +AD
					FREQ = FreqD[ALT]

					sampSiteInfoL = []
					sampSiteInfoL.append(REF)
					sampSiteInfoL.append(ALT)
					sampSiteInfoL.append(MDP)
					sampSiteInfoL.append(DP)
					sampSiteInfoL.append(RD)
					sampSiteInfoL.append(AD)
					sampSiteInfoL.append(FREQ)
					sampSiteInfo = "\t".join(sampSiteInfoL)
					sampSiteInfoD[Site] = sampSiteInfo

				else:
					notHegeL.append(Site)
	return sampSiteInfoD				



def vcfFomat(sampPosL,snvInfoD,SNPInfoD,RefChromID):
	outL = []
	NotHegeCount = 0
	for sampPos in sampPosL:
		lineL = ''
		if sampPos in snvInfoD.keys():
			lineL = snvInfoD[sampPos].split("\t")
		elif sampPos in SNPInfoD.keys():
			lineL = SNPInfoD[sampPos].split("\t")
		else:
			NotHegeCount += 1

		
		outlineL = []
		if lineL != '':
			outlineL.append(RefChromID)
			outlineL.append(sampPos)
			outlineL.append(".")
			outlineL.append(lineL[0])  ##REF
			outlineL.append(lineL[1])  ##ALT
			outlineL.append(".")
			outlineL.append("PASS")
			outlineL.append(";".join(["ADP="+lineL[3],"WT=0","HET=0","HOM=1","NC=0"]))
			outlineL.append("GT:GQ:MDP:DP:RD:AD:FREQ")
			outlineL.append(":".join(["1/1","255",lineL[2],lineL[3],lineL[4],lineL[5],lineL[6]+"%"]))
			
			outline = "\t".join(outlineL)
			outL.append(outline)
	out = "\n".join(outL)
	
	return out




