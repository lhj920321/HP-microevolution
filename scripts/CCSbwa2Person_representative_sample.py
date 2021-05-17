###################################
###################################
#######CCS data map to respe genome
###################################

#person = "P2 P18 P25 P5"

person = "P25"


if person == 'P1':
	bwa_ref_samp = 'P1-E-j'
	REP_INDEX=["P1-E-x","P1-E-d","P1-E-j","P1-E-t"]

if person == 'P2':
	bwa_ref_samp = "P2-E-t"
	REP_INDEX=["P2-E-t"]

if person == 'P4':
	bwa_ref_samp = "P4-E-j"
	REP_INDEX= ["P4-E-x","P4-E-d","P4-E-j"]	

if person == 'P5':
	bwa_ref_samp =  "P5-E-d"
	REP_INDEX= ["P5-E-d"]	

if person == 'P6':
	bwa_ref_samp = "P6-E-t"
	REP_INDEX= ["P6-E-d","P6-E-t"]	

if person == 'P7':
	bwa_ref_samp = "P7-E-j"
	REP_INDEX= ["P7-E-d","P7-E-j"]	


if person == 'P18':
	bwa_ref_samp = "P18-C_bR-x"
	REP_INDEX= ["P18-C_bR-x","P18-C_bR-j"]	

if person == 'P19':
	bwa_ref_samp = "P19-C_bR-j"
	REP_INDEX= ["P19-C_bR-d","P19-C_bR-j","P19-C_bR-t"]	

if person == 'P20':
	bwa_ref_samp = "P20-C_bR-t" 
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
	## bwa to self  genome
		expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.bam",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),
		expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.bam.bai",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp),

		#expand("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/varscan_snv/{sample}_2_{bwa_ref_samp}.sort.varscan2.mpileup2snp.FreqStat.pdf",sample=REP_INDEX,bwa_ref_samp=bwa_ref_samp)







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
		"/shared/liuhj/HP/data/CCS/{sample}.fastq.gz",
		"/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta.bwt"
	output:
		temp("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/bwa/{sample}/{sample}_2_{bwa_ref_samp}.bam")
	shell:
		"bwa  mem  -t  30  -x  pacbio  {input[0]}  {input[1]}  \
      		|  samtools view -bS -bF 4 -  > {output}   "   #
 

rule bam_sort:
	input:
		temp("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/bwa/{sample}/{sample}_2_{bwa_ref_samp}.bam"),
	output:
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.bam",
	shell:
		"samtools  sort  --threads   30  {input}  -o {output}"





rule bam_index:
	input:
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.bam"
	output:
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/bwa/{sample}/{sample}_2_{bwa_ref_samp}.sort.bam.bai"
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
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/bwa/{sample}/{sample}_2_{bwa_ref_samp}_2_selfGnm.sort.bam",
		"/shared/liuhj/HP/process/assembly/all_genomes/{bwa_ref_samp}.fasta.fai",
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/bwa/{sample}/{sample}_2_{bwa_ref_samp}_2_selfGnm.sort.bam.bai"
	output:
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/mpileup/all_mpileup/{sample}_2_{bwa_ref_samp}.sort.mpileup"
	log:
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/mpileup/all_mpileup/{sample}_2_{bwa_ref_samp}.sort.mpileup.log"
	shell:
		"samtools mpileup -q 30  --reference  {input[0]}  {input[1]}  1>{output} 2>{log}"     

#-Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
#-q, --min-MQ INT        skip alignments with mapQ smaller than INT [0]
#-d, –max-depth 最大测序深度，过滤掉超深度测序的位点




rule varscan2_mpileup2snp:
	input:
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/mpileup/all_mpileup/{sample}_2_{bwa_ref_samp}.sort.mpileup"
	output:
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/varscan_snv/{sample}_2_{bwa_ref_samp}.sort.varscan2.mpileup2snp.vcf"
	shell:
		"varscan   mpileup2snp  {input}  --min-coverage 50  --min-reads2  5  --min-avg-qual 30  --min-var-freq  0.02   --output-vcf 1  >{output} "     


rule varscan_snv_freq_distribution:
	input:	
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/varscan_snv/{sample}_2_{bwa_ref_samp}.sort.varscan2.mpileup2snp.vcf",
	output:
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/varscan_snv/{sample}_2_{bwa_ref_samp}.sort.varscan2.mpileup2snp.FreqStat.txt",
	params:
		minFreq="0",
		maxFreq="100",
		frameStep="2",
		scriptP="/shared/liuhj/HP/scripts/pipeline_scripts/iSNVcalling_snpeff_FreqStat"
	shell:
		"python {params[3]}/varscan_freq_dictribution.py  -i {input}  -o  {output} -m  {params[0]}  -M  {params[1]}  -s   {params[2]}  "  



rule varscan_snv_freq_distribution_plot:
	input:	
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/varscan_snv/{sample}_2_{bwa_ref_samp}.sort.varscan2.mpileup2snp.FreqStat.txt",
	output:
		"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/varscan_snv/{sample}_2_{bwa_ref_samp}.sort.varscan2.mpileup2snp.FreqStat.pdf",
	params:
		scriptP="/shared/liuhj/HP/scripts/pipeline_scripts/iSNVcalling_snpeff_FreqStat"
	shell:
		"python {params[0]}/varscan_freq_dictribution_barplot.py  -i {input}  -o  {output}"  




