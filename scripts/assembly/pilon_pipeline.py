###################################
#######pilon pipeline 
#######
###################################

REP_INDEX = ["P12-E-x"]


rule all:
	input:
	##trimmmomatic
		#expand("/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R1.paired.fastq.gz",sample=REP_INDEX),
		#expand("/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R2.paired.fastq.gz",sample=REP_INDEX),
	##pilon 
		expand("/shared/liuhj/HP/process/assembly/ONT_assembly/pilon/{sample}/{sample}.pilon.1.fasta",sample=REP_INDEX)


rule reads_Trim:
	input:
		"/shared/liuhj/HP/data/NGS/LBFC20201376/210109_A00838_0374_BHV3CTDSXY/{sample}.R1.fastq.gz",
		"/shared/liuhj/HP/data/NGS/LBFC20201376/210109_A00838_0374_BHV3CTDSXY/{sample}.R2.fastq.gz",
	params:
		"/shared/liuhj/HP/data/NGS/LBFC20201376/trimmomatic",
	threads:8
	output:
		"/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R1.paired.fastq.gz",
		temp("/shared/liuhj/HP/data/NGS/LBFC20201376/trimmomatic/{sample}.R1.unpaired.fastq.gz"),
		"/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R2.paired.fastq.gz",
		temp("/shared/liuhj/HP/data/NGS/LBFC20201376/trimmomatic/{sample}.R2.unpaired.fastq.gz"),
	log:
		temp("/shared/liuhj/HP/data/NGS/LBFC20201376/trimmomatic/{sample}.trim.log"),
		temp("/shared/liuhj/HP/data/NGS/LBFC20201376/trimmomatic/{sample}.trim.summary.txt"),
	shell:
		"trimmomatic  PE   -threads  30   {input[0]} {input[1]} \
		{output[0]} {output[1]} {output[2]} {output[3]}  -trimlog  {log[0]}  -summary  {log[1]} \
 		SLIDINGWINDOW:4:20  HEADCROP:10 "   



rule pilon:
	input:	
		"/shared/liuhj/HP/process/assembly/ONT_assembly/HP_P12-E-x/6-bridge_contigs/bridged_contigs.fasta",   ##需要改
		"/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R1.paired.fastq.gz",  ##需要改
		"/shared/liuhj/HP/data/NGS/Trimmomatic/{sample}/{sample}.R2.paired.fastq.gz",  ##需要改
	output:
		"/shared/liuhj/HP/process/assembly/ONT_assembly/pilon/{sample}/{sample}.pilon.1.fasta",
		temp("/shared/liuhj/HP/process/assembly/ONT_assembly/pilon/{sample}/{sample}.0.bam"),
		temp("/shared/liuhj/HP/process/assembly/ONT_assembly/pilon/{sample}/{sample}.1.bam"),
		temp("/shared/liuhj/HP/process/assembly/ONT_assembly/pilon/{sample}/{sample}.2.bam"),
	params:
		sample="{sample}",
		step="pilon",
		outputP="/shared/liuhj/HP/process/assembly/ONT_assembly/pilon/{sample}",
		pilon_times=1,
	log:
		"/shared/liuhj/HP/process/assembly/ONT_assembly/pilon/{sample}/{sample}.pilon.log",
	shell:
		"bash  Pilon.sh   {params[0]}  {params[1]} \
		{input[0]}   {input[1]}  {input[2]}  {params[2]}   {params[3]}   "  








