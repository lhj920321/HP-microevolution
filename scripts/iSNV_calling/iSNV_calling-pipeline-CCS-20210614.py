#coding=UTF-8
###################################
#######iSNV_calling 20210517
#######
###################################

#CCS person
personLst = ["P1","P2","P4","P5","P6","P7","P18","P19","P20","P21","P22","P23","P24","P25"]

#personLst = ["P1"]


'''

STAR_INDEX = "/app/ref/ensembl/human/idx/" 
SAMPLE_DIR = "Rawdata" 
SAMPLES, = glob_wildcards(SAMPLE_DIR + "/{sample}_1.fq.gz") 

R1 = '{sample}_1.fq.gz' 
R2 = '{sample}_2.fq.gz' 

rule all: 
    input: 
     expand("Align/{sample}/Aligned.toTranscriptome.out.bam", sample=SAMPLES) 

rule alignment: 
    input: 
     r1 = join(SAMPLE_DIR, R1), 
     r2 = join(SAMPLE_DIR, R2) 
    params: 
     STAR_INDEX = STAR_INDEX 
    output: 
     "Align/{sample}/Aligned.toTranscriptome.out.bam" 
    threads: 
     8 
    message: 
     "--- Mapping STAR---" 
    shell:""" 
'''



rule all:
	input:
	#bwa map
		## iSNV calling
		expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/all.iSNV_with_SNP.pyResults.txt",person=personLst),
		##vcf
		expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/vcf",person=personLst),
		##snpeff vcf
		#expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/vcf_snpeff",person=personLst),

		##general snvFreq_distribut_Stat
		#expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/general-snvFreq_distribut_Stat",person=personLst),
		## snvNum_Stat
		#expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/table/snvNum_Stat",person= personLst),
		## snv_density_distribution
		#expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/snv_density_Stat_1k",person= personLst),
		
		##all sa -----
		#expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/20210331-snvSNP_AnnoStat_OtherBacFilted_RepeatFilted",person=personLst),

		##SNP genome sub SNP posi
		#expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/subGnm",person=personLst),
		#expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/subGnm_prokka_gff",person=personLst),


		##SNV Freq distribut(Repeat Regions were filted)
		expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/SNVFreq_distribut_Stat_filtRepeat",person=personLst),




rule ntfreq_2_FreqBigtable:
	input:
		configF="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}.ntfreq_file_list.txt",  ##format in the config file: path/sampID.ntfreq"\n"
		refF="/shared/liuhj/HP/process/assembly/person_respectGnm/{person}.respect.fasta",
	output:	
		outputF="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/all.iSNV_with_SNP.pyResults.txt",
		outputF2="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/out.samps-ntfreq.bigtable.txt",

	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		outputP="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}",
		DEP_THRES=100,       ##最低深度
		MIN_DEP_THRES=5,    ##支持alt的最小reads数
		VALID_SIZE_THRES=1000000,  ##样本有效需要的有效位点数
		FREQ_THRES=0.02,    ###alt的频率最小值
		STRANDED_RATIO_THRES=0.1,     ## reads链偏性检测

	shell:
		"python {params.scriptPath}/step2_ntfreq_To_bigtables.py  -i  {input.configF} -r {input.refF}  -o {params.outputP}  \
		-D  {params.DEP_THRES} -A {params.MIN_DEP_THRES} -N {params.VALID_SIZE_THRES} -F {params.FREQ_THRES} -S {params.STRANDED_RATIO_THRES} "   #




rule iSNV_freq_distribut:
	input:
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/all.iSNV_with_SNP.pyResults.txt",
	output:
		outputP=directory("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/general-snvFreq_distribut_Stat")
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		#iSNV_FreqThres = 0.05,     ## for examp: iSNV_FreqThres = 0.05, mean iSNV freq :>=0.05 && <0.95  ;SNP: >=0.95
		minFreq	= 0.02,
		maxFreq	= 1.0,
		step = 2,

	shell:
		"python  {params.scriptPath}/iSNVpy_Freq_distribut_Stat.py  -i {input} -o  {output.outputP} -m {params.minFreq}  -M {params.maxFreq}  -s {params.step}  "





rule iSNV_freq_distribut_RepeatRegionfilt:
	input:
		iSNVtable="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/all.iSNV_with_SNP.pyResults.txt",
		#ClonalFrameF="/shared/liuhj/HP/process/samps_Tree-20210527/clonalFrame.RAxML_bestTree.20210527-allSamp.raxml.importation_status.txt",
		RepeatRegionP="/shared/liuhj/HP/process/Genome_repeatRegion",
	output:
		outputP=directory("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/SNVFreq_distribut_Stat_filtRepeat")
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		#iSNV_FreqThres = 0.05,     ## for examp: iSNV_FreqThres = 0.05, mean iSNV freq :>=0.05 && <0.95  ;SNP: >=0.95
		minFreq	= 0.05,
		maxFreq	= 1.0,
		step = 2,
		seqTech = "CCS",

	shell:
		"python  {params.scriptPath}/iSNVpy_Freq_distribut_Stat-filtRepeatRegion.py  -i {input} -o  {output.outputP}  -R {input.RepeatRegionP} -m {params.minFreq}  -M {params.maxFreq}  -s {params.step} -T {params.seqTech} "






rule convert_iSNVpytable_To_vcf :
	input:
		iSNVtable="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/all.iSNV_with_SNP.pyResults.txt",
		ntfreqTable = "/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/out.samps-ntfreq.bigtable.txt",
		RefGnmF = "/shared/liuhj/HP/process/assembly/person_respectGnm/{person}.respect.fasta"

	output:
		outputP=directory("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/vcf")
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		MinFreq	= 0.02,
		#snv_MaxFreq	= 0.95,
		

	shell:
		"python  {params.scriptPath}/iSNVpy_iSNVtable_To_vcf.py   -i {input.iSNVtable}  \
		-n {input.ntfreqTable}  -r  {input.RefGnmF} -o  {output.outputP} -m {params.MinFreq}  "


##snpeff ann

#for file in ../../process/person_ntfreq_tables/P*/vcf/*.vcf ;do echo $file ;filename=${file##*/};fileID=${filename%.vcf*}; snpEff  ann   -v  AE004091  $file > ${file%/*}/$fileID.snpeff.vcf ;done

rule snpeffAnno:
	input:
		inputvcfP="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/vcf"
	output:
		outputP=directory("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/vcf_snpeff")
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		dbID="AE004091",
		snpeffSoftP="/shared/liuhj/software/snpEff",
	shell:
		"bash {params.scriptPath}/pipeline-snpeff.sh  {input} {output} {params.dbID} {params.snpeffSoftP}"

'''
vcfP=$1
outSnpeffVcfP=$2
dbID=$3
snpeffSoftP=$4
'''






rule iSNV_snvNumCount:
	input:
		iSNVtable="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/all.iSNV_with_SNP.pyResults.txt",
		infoF = "/shared/liuhj/tonglv/Infos/20210311-info.txt",
	output:
		outputP=directory("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/snvNum_Stat")
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		#iSNV_FreqThres = 0.05,     ## for examp: iSNV_FreqThres = 0.05, mean iSNV freq :>=0.05 && <0.95  ;SNP: >=0.95
		minFreq	= 0.1,
		maxFreq	= 0.9,
		

	shell:
		"python  {params.scriptPath}/iSNVpy_Freq_distribut_Stat.py  -i {input.iSNVtable}  \
		-f {input.infoF}  -o  {output.outputP} -m {params.minFreq}  -M {params.maxFreq}  "


rule iSNV_snv_density_Stat_1k:
	input:
		iSNVtable="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/all.iSNV_with_SNP.pyResults.txt",
	output:
		outputP=directory("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/snv_density_Stat_1k"),
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		minFreq	= 0.02,
		#minFreq for SNV ,if = 0.05,     ## for examp: iSNV_FreqThres = 0.05, mean iSNV freq :>=0.05 && <0.95  ;SNP: >=0.95
		denst_step = 500,
		BacRefLen = 1700000,
	shell:
		"python iSNVpy_SNV_density_distribut.py -i  {input.iSNVtable} -o  {output.outputP} -M {params.minFreq} -S {params.denst_step} -L {params.BacRefLen}"


'''
rule iSNV_snvNumCount_noSNV:
	input:
		iSNVtable="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/all.iSNV_with_SNP.pyResults.txt",
		infoF = "/shared/liuhj/tonglv/Infos/20210311-info.txt",
		filtHomoPosi="/shared/liuhj/tonglv/process/20210402-bwaToOtherBacteria/PA-HomoRegions-FromOtherBacteria.txt",
		snpEffPath="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/vcf_snpeff",
		RepeatMaskedGnm="/shared/liuhj/tonglv/process/RefGnmRepeatRegion/repeatmasker_out/AE004091.2.fasta.masked.fasta",
		RepeatPosiF="/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/PA.RepeatRegion-Cov30-identy60.Posi.txt"
	output:
		outputP=directory("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/table/snvNum_Stat")
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		#iSNV_FreqThres = 0.05,     ## for examp: iSNV_FreqThres = 0.05, mean iSNV freq :>=0.05 && <0.95  ;SNP: >=0.95
		minFreq	= 0.1,
		maxFreq	= 0.9,
		
	shell:
		"python  {params.scriptPath}/iSNVpy_snvSNP_Stat_clade.py   -i {input.iSNVtable}  \
		-f {input.infoF} -e {input.snpEffPath} -R {input.RepeatPosiF} -o  {output.outputP}  -m {params.minFreq}  -M {params.maxFreq}  -E {input.filtHomoPosi}   "

'''

rule iSNV_SNP_consensusGnm:
	input:
		iSNVtable="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/all.iSNV_with_SNP.pyResults.txt",
		ntfreqTable = "/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/out.samps-ntfreq.bigtable.txt",
		RefGnmF = "/shared/liuhj/HP/process/assembly/person_respectGnm/{person}.respect.fasta"

	output:
		outputP=directory("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/subGnm")
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		MinFreq	= 0.9,
		#snv_MaxFreq	= 0.95,
		

	shell:
		"python  {params.scriptPath}/iSNVpy_iSNVtable_RefGnmtihuanSNP.py   -i {input.iSNVtable}  \
		-n {input.ntfreqTable}  -r  {input.RefGnmF} -o  {output.outputP} -m {params.MinFreq}  "






rule consensusGnm_prokkaAnno:
	input:
		inputP="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/subGnm"
	output:
		outputP=directory("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/{person}/subGnm_prokka_gff")
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		genbankdb="/home/amax/anaconda3/envs/biotools_lhj/db/Pseudomonas/AE004091.2.gb"

	shell:
		"{params.scriptPath}/pipeline-prokka.sh  {input.inputP}  {output.outputP}  {params.genbankdb} "


#gffP=$1
#outP=$2
#genbankdb=$3  ## tonglv --- /home/amax/anaconda3/envs/biotools_lhj/db/Pseudomonas/AE004091.2.gb











