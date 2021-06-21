###################################
###################################
#######
###################################

#person_list="7 11 12    21 201 213 13    124 1346 142 75    340 345 144 457"

person = "P11"


if person == 'P1':
	bwa_ref_samp = 'P1-E-j'
	REP_INDEX=["P1-E-x","P1-E-d","P1-E-j","P1-E-t"]

if person == 'P2':
	bwa_ref_samp = "P2-E-t"
	REP_INDEX=["P2-E-x","P2-E-d","P2-E-j","P2-E-t"]

if person == 'P3':
	bwa_ref_samp = "P3-E-j"
	REP_INDEX= ["P3-E-x","P3-E-j","P3-E-t"]

if person == 'P4':
	bwa_ref_samp = "P4-E-j"
	REP_INDEX= ["P4-E-x","P4-E-d","P4-E-j"]

if person == 'P5':
	bwa_ref_samp = "P5-E-d"
	REP_INDEX= ["P5-E-d","P5-E-j"]

if person == 'P6':
	bwa_ref_samp = "P6-E-t"
	REP_INDEX= ["P6-E-d","P6-E-t"]

if person == 'P7':
	bwa_ref_samp = "P7-E-j"
	REP_INDEX= ["P7-E-d","P7-E-j","P7-E-t"]

if person == 'P8':
	bwa_ref_samp = "P8-E-x"
	REP_INDEX= ["P8-E-t"]

if person == 'P9':
	bwa_ref_samp = "P9-E-j"
	REP_INDEX= ["P9-E-t"]

if person == 'P10':
	bwa_ref_samp = "P10-E-j"
	REP_INDEX= ["P10-E-t"]

'''
if person == 'P8':
	bwa_ref_samp = "P8-E-x"
	REP_INDEX= ["P8-E-x","P8-E-d","P8-E-t"]

if person == 'P9':
	bwa_ref_samp = "P9-E-j"
	REP_INDEX= ["P9-E-d","P9-E-j","P9-E-t"]

if person == 'P10':
	bwa_ref_samp = "P10-E-j"
	REP_INDEX= ["P10-E-x","P10-E-j","P10-E-t"]
'''

if person == 'P11':
	bwa_ref_samp = "P11-E-x"
	#REP_INDEX= ["P11-E-x","P11-E-j","P11-E-t"]
	REP_INDEX= ["P11-E-x","P11-E-t"]


if person == 'P12':
	bwa_ref_samp = "P12-E-t"
	REP_INDEX= ["P12-E-x","P12-E-t"]

if person == 'P13':
	bwa_ref_samp = "P13-E-j"
	REP_INDEX= ["P13-E-x","P13-E-d","P13-E-j","P13-E-t"]

if person == 'P14':
	bwa_ref_samp = "P14-E-j"
	REP_INDEX= ["P14-E-x","P14-E-d","P14-E-j"]

if person == 'P15':
	bwa_ref_samp = "P15-C_bS-j"
	REP_INDEX= ["P15-C_bS-d","P15-C_bS-j","P15-C_bS-t"]

if person == 'P16':
	bwa_ref_samp = "P16-C_bS-j"
	REP_INDEX= ["P16-C_bS-x","P16-C_bS-d","P16-C_bS-j"]

if person == 'P17':
	bwa_ref_samp = "P17-C_bR-j"
	REP_INDEX= ["P17-C_bR-x","P17-C_bR-d","P17-C_bR-j"]

if person == 'P18':
	bwa_ref_samp = "P18-C_bR-x"
	REP_INDEX= ["P18-C_bR-x","P18-C_bR-j"]

if person == 'P19':
	bwa_ref_samp = "P19-C_bR-j"
	REP_INDEX= ["P19-C_bR-d","P19-C_bR-j","P19-C_bR-t"]

if person == 'P20':
	bwa_ref_samp = "P20-C_bR-t"  ##changed
	REP_INDEX= ["P20-C_bR-x","P20-C_bR-d","P20-C_bR-j","P20-C_bR-t"]

if person == 'P21':
	bwa_ref_samp = "P21-C_sL-j"
	REP_INDEX= ["P21-C_sL-x","P21-C_sL-d","P21-C_sL-j","P21-C_sL-t"]

if person == 'P22':
	bwa_ref_samp = "P22-C_sL-d"  ##changed
	REP_INDEX= ["P22-C_sL-d","P22-C_sL-t"]

if person == 'P23':
	bwa_ref_samp = "P23-C_sL-j"
	REP_INDEX= ["P23-C_sL-x","P23-C_sL-j"]

if person == 'P24':
	bwa_ref_samp = "P24-C_sC-t"  ##changed
	REP_INDEX= ["P24-C_sC-x","P24-C_sC-j","P24-C_sC-t"]

if person == 'P25':
	bwa_ref_samp = "P25-C_sC-t"
	REP_INDEX= ["P25-C_sC-j","P25-C_sC-t"]




rule all:
	input:
	## trimmomatic   这一步就别关闭了，后边有的程序要用
		expand("/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R1.paired.fastq.gz",sample=REP_INDEX),
		expand("/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R1.paired.fastq.gz",sample=REP_INDEX),
	## NECAT  (在smp2上)
		#expand("process/NECAT/{sample}/6/{sample}.bridged_contigs.fasta",sample=REP_INDEX),
	## pilon   --1
		#expand("process/pilon/{bwa_ref_samp}/{bwa_ref_samp}.pilon.3.fasta",sample=REP_INDEX),
	##circlator
		#expand("process/circlator/{sample}/04.merge.circularise_details.log",sample=REP_INDEX),
	## bwa to representive  genome
		expand("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.bam",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),
		expand("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.bam.bai",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),

	## bam to mpileup
		expand("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/mpileup/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.mpileup",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),
	## varscan2
		expand("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/mpileup/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.varscan2.mpileup2snp.vcf",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),
	##varscan_Freq_distrib
		#expand("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/varscan2_FreqDistribut/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.varscan2.mpileup2snp.FreqStat.txt",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),
	##varscan_Freq_distrib  bar plot
		#expand("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/varscan2_FreqDistribut/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.varscan2.mpileup2snp.FreqStat.pdf",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),

	## pindel
		#expand("process/pindel/BGI_2_PersonReprest_genome/{sample}/{sample}_2_{bwa_ref_samp}.confg.txt",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),
		#expand("process/pindel/BGI_2_PersonReprest_genome/{sample}/{sample}_2_{bwa_ref_samp}.pindel.Indel.vcf",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),
	##indel filter
		#expand("process/common_snv_indel/BGI_2_PersonReprest_genome/filted_indel/{sample}_2_{bwa_ref_samp}.pindel.Indel_filted.vcf",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),


	## prokka  基因组功能注释
		#expand("process/prokka/anno_used_26695/{sample}/{sample}.gff",sample=REP_INDEX),
		#expand("process/prokka/anno_used_26695/{sample}/{sample}.ffn",sample=REP_INDEX),
		#expand("process/prokka/anno_used_26695/{sample}/{sample}.tsv",sample=REP_INDEX),
		#expand("process/prokka/anno_used_26695/{sample}/{sample}.faa",sample=REP_INDEX),
	##snpEff --snv (在windows系统的虚拟机中)
		##运行 snpEff_snv.sh	
	##snpEff --indel (在windows系统的虚拟机中)	
		##运行 snpEff_indel.sh			

	##ngmlr  ONT比对到代表基因组上
		#expand("process/nglmr/ONT_2_PersonReprest_genome/{sample}/ONT_{sample}_nglmr2_{bwa_ref_samp}.bam",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),
		#expand("process/nglmr/ONT_2_PersonReprest_genome/{sample}/ONT_{sample}_nglmr2_{bwa_ref_samp}.sort.Header.bam",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),
		#expand("process/nglmr/ONT_2_PersonReprest_genome/{sample}/ONT_{sample}_nglmr2_{bwa_ref_samp}.sort.Header.bam.bai",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),

	##sniffles  专门找三代数据的SV
		#expand("process/sniffles/ONT_2_PersonReprest_genome/{sample}/ONT_{sample}_nglmr2_{bwa_ref_samp}.sort.bam.sniffles.vcf",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),

##-----统计分析
	## 统计snv的freq 分布，看是否有峰存在
		##运行 common_snv_freq_dictribution.py  $person
		##运行 R脚本画图  snv_freq_hist.r


'''
rule reads_Trim:
	input:
		"/shared/liuhj/HP/data/NGS/LBFC20201376/210109_A00838_0374_BHV3CTDSXY/{sample}.NGS.R1.fastq",
		"/shared/liuhj/HP/data/NGS/LBFC20201376/210109_A00838_0374_BHV3CTDSXY/{sample}.NGS.R2.fastq",
	params:
		"/shared/liuhj/HP/data/NGS/trimmomatic",
	threads:8
	output:
		"/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R1.paired.fastq.gz",
		"/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}.R1.unpaired.fastq.gz",
		"/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R2.paired.fastq.gz",
		"/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}.R2.unpaired.fastq.gz",
	log:
		temp("/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}.trim.log"),
		temp("/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}.trim.summary.txt"),
	shell:
		"trimmomatic  PE   -threads  30   {input[0]} {input[1]} \
		{output[0]} {output[1]} {output[2]} {output[3]}  -trimlog  {log[0]}  -summary  {log[1]} \
 		SLIDINGWINDOW:4:20  HEADCROP:10 TRAILING:30"   





rule circlator_cir:
	input:
		"process/pilon/{bwa_ref_samp}/{bwa_ref_samp}.pilon.3.fasta",
		"data/ONT/barcode/{sample}.fastq"	  #"data/ONT/barcode/third/{sample}.fastq"	
	output:
		"process/circlator/{sample}/04.merge.circularise_details.log"
	params:
		"process/circlator/{sample}",
		#merge_min_id=85,
		#merge_breaklen=1000,
	shell:
		"rm -r {params[0]} "
		"circlator all --merge_min_id  85 --merge_breaklen 1000  {input[0]}  {input[1]}  {params[0]}"



'''


rule bwa_index:
	input:
		"/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta"
	output:
		"/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta.bwt"
	shell:
		"bwa index {input}"


rule bwa_run:
	input:
		"/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta",
		"/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R1.paired.fastq.gz",
		"/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R2.paired.fastq.gz",
		"/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta.bwt"
	output:
		temp("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.bam")
	shell:
		"bwa  mem  -t  30  {input[0]}  {input[1]}  {input[2]} \
      		|  samtools view -bS -  > {output}   "   #
 

rule bam_sort:
	input:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.bam"
	output:
		temp("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.bam"),
	shell:
		"samtools  sort  --threads   30  {input}  -o {output}"




      	
rule bam_rmdup:
	input:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.bam"
	output:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.bam",
		temp("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.metrics")
	log:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.log"
	shell:
		"picard  MarkDuplicates  INPUT={input}  \
		OUTPUT={output[0]}  METRICS_FILE={output[1]}   REMOVE_DUPLICATES=true  2> {log}"



rule BGI_bam_Header:
	input:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.bam"
	output:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.bam"
	params:
		ID_SM="sample_{sample}",
		LB_name="hp_{sample}_lib1"
	shell:
		"picard  AddOrReplaceReadGroups I={input}  O={output} \
		ID={params[0]} LB={params[1]} PL=NGS  SM={params[0]}   PU=unit1"



rule bam_index:
	input:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.bam"
	output:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.bam.bai"
	shell:
		"samtools index  {input}"




rule refGonm_faidx:
	input:
		"/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta"
	output:
		"/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta.fai"
	shell:
		"samtools  faidx  {input}"




rule mpileup:
	input:
		"/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta",
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.bam",
		"/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta.fai",
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.bam.bai"
	output:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/mpileup/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.mpileup"
	log:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/mpileup/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.mpileup.log"
	shell:
		"samtools mpileup -q 30  --reference  {input[0]}  {input[1]}  1>{output} 2>{log}"     

#-Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
#-q, --min-MQ INT        skip alignments with mapQ smaller than INT [0]
#-d, –max-depth 最大测序深度，过滤掉超深度测序的位点


rule varscan2_mpileup2snp:
	input:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/mpileup/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.mpileup"
	output:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/mpileup/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.varscan2.mpileup2snp.vcf"
	shell:
		"varscan   mpileup2snp  {input}  --min-coverage 50  --min-reads2  5  --min-avg-qual 30  --min-var-freq  0.02   --output-vcf 1  >{output} "     


rule varscan_snv_freq_distribution:
	input:	
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/mpileup/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.varscan2.mpileup2snp.vcf",
	output:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/varscan2_FreqDistribut/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.varscan2.mpileup2snp.FreqStat.txt",
	params:
		minFreq="0",
		maxFreq="100",
		frameStep="2",
		scriptP="/shared/liuhj/HP/scripts/pipeline_scripts/iSNVcalling_snpeff_FreqStat"
	shell:
		"python {params[3]}/varscan_freq_dictribution.py  -i {input}  -o  {output} -m  {params[0]}  -M  {params[1]}  -s   {params[2]}  "  


rule varscan_snv_freq_distribution_plot:
	input:	
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/varscan2_FreqDistribut/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.varscan2.mpileup2snp.FreqStat.txt",
	output:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/varscan2_FreqDistribut/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.varscan2.mpileup2snp.FreqStat.pdf",
	params:
		scriptP="/shared/liuhj/HP/scripts/pipeline_scripts/iSNVcalling_snpeff_FreqStat"
	shell:
		"python {params[0]}/varscan_freq_dictribution_barplot.py  -i {input}  -o  {output}"  



rule pindel_indel:
	input:
		ref_F="/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta",
		bam_F="/shared/liuhj/HP/process/NGS_map2PersonResptGnm/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.bam",
	output:	
		"process/pindel/BGI_2_PersonReprest_genome/{sample}/{sample}_2_{bwa_ref_samp}.confg.txt",  ##pindel_confg_F=
		"process/pindel/BGI_2_PersonReprest_genome/{sample}/{sample}_2_{bwa_ref_samp}.pindel.Indel.vcf",  #pindel_indel_F=
	params:
		"{sample}",
		chrom_NAME="{bwa_ref_samp}_pilon_pilon_pilon",
		pindel_outdir="process/pindel/BGI_2_PersonReprest_genome/{sample}",
		min_mapping_Q=30,
		min_alt_reads=10
	log:
		"process/pindel/BGI_2_PersonReprest_genome/{sample}/{sample}_2_{bwa_ref_samp}.pindel.Indel.log"
	shell:
		"bash pindel.sh  {input[0]} {input[1]} {output[0]} {output[1]} {params[0]} {params[1]} {params[2]} {params[3]} {params[4]}  2>{log} "

#ref_F=$1
#bam_F=$2
#pindel_confg_F=$3
#pindel_indel_F=$4
#sample=$5
#chrom_NAME=$6
#pindel_outdir=$7
#minimal_mapping_quality=30
#min_alt_reads=$9
# -T/--number_of_threads the number of threads Pindel will use (default 1).


rule pindel_filter:
	input:
		pindel_snp_vcf="process/pindel/BGI_2_PersonReprest_genome/{sample}/{sample}_2_{bwa_ref_samp}.pindel.Indel.vcf",
	output:	
		out_pindel_indel="process/common_snv_indel/BGI_2_PersonReprest_genome/filted_indel/{sample}_2_{bwa_ref_samp}.pindel.Indel_filted.vcf",  ##pindel_confg_F=
	params:
		"{sample}",
		"{bwa_ref_samp}",
		max_same_char_num=3,  ##控制位点前后的重复碱基和重复二碱基的数量不超过多少
		min_indel_freq=10,  ##10%
		min_Posireads=80,  ##该位点的最小reads数

	shell:
		"python indel_filter.py  {params[0]} {params[1]} {params[2]} {params[3]} {params[4]}  {input[0]} {output[0]} "

#sample = str(sys.argv[1])
#ref_samp = str(sys.argv[2])
#max_same_char_num=int(sys.argv[3]),
#min_indel_freq=int(sys.argv[4]),  ##10%
#min_Posireads=int(sys.argv[5]),
#pindel_snp_vcf = str(sys.argv[6])
#out_pindel_indel = str(sys.argv[7]) 





rule varscanIndel:
	input:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/mpileup/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.mpileup"
	output:
		"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/mpileup/{sample}/{sample}_2_{bwa_ref_samp}.sort.rmdup.Header.mpileup2Indel.vcf"
	shell:
		"varscan   mpileup2indel  {input}  --min-coverage 100  --min-reads2  10  --min-avg-qual 30  --min-var-freq  0.1   --output-vcf 1  >{output}"




