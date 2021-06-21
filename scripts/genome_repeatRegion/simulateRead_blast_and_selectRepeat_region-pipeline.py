###################################
#######1.
#######
###################################

Samps=["P10-E-j","P10-E-t","P10-E-x","P11-E-j.ONT","P11-E-t","P11-E-x","P12-E-t","P12-E-x","P13-E-d","P13-E-j","P13-E-x","P14-E-d","P14-E-j","P14-E-x","P15-C_bS-d","P15-C_bS-j","P15-C_bS-t","P16-C_bS-d","P16-C_bS-j","P16-C_bS-x","P17-C_bR-d","P17-C_bR-j","P17-C_bR-x","P18-C_bR-j","P18-C_bR-x","P19-C_bR-d","P19-C_bR-j","P19-C_bR-t","P1-E-d","P1-E-j","P1-E-t","P1-E-x","P20-C_bR-d","P20-C_bR-j","P20-C_bR-t","P20-C_bR-x","P21-C_sL-d","P21-C_sL-j","P21-C_sL-t","P21-C_sL-x","P22-C_sL-d","P22-C_sL-t","P23-C_sL-j","P23-C_sL-x","P24-C_sC-j","P24-C_sC-t","P24-C_sC-x","P25-C_sC-j","P25-C_sC-t","P2-E-d","P2-E-j","P2-E-t","P2-E-x","P3-E-j","P3-E-t","P3-E-x","P4-E-d","P4-E-j","P4-E-x","P5-E-d","P5-E-j","P6-E-d","P6-E-t","P7-E-d","P7-E-j","P7-E-t","P8-E-d","P8-E-t","P8-E-x","P9-E-d","P9-E-j","P9-E-t"]

#Samps=["P10-E-t"]



rule all:
	input:
		expand("/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.simulateReads.fasta",samp=Samps),

		expand("/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.RepeatRegion.txt",samp=Samps),

		#expand("/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.blastOut.txt",samp=Samps),



rule cutGenmTosimulRead:
	input:
		GnmF="/shared/liuhj/HP/process/assembly/all_genomes/{samp}.fasta"
	output:
		GnmSimulReadsF="/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.simulateReads.fasta"
	params:
		scriptP="/shared/liuhj/HP/scripts/HP-microevolution/scripts/genome_repeatRegion",
		cutStep="75",
		cutlen="150",	
	shell:
		"python  {params.scriptP}/RefGnm_cutToSimmulateReads.py  -r  {input.GnmF}  -o {output.GnmSimulReadsF} -s {params.cutStep} -l {params.cutlen} "



rule makeblastdb:
	input:
		GnmF="/shared/liuhj/HP/process/assembly/all_genomes/{samp}.fasta"
	output:
		GnmBlasdbF="/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.blastdb.ndb"
	params:
		GnmBlasdbID="/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.blastdb"


	shell:
		"makeblastdb  -in   {input.GnmF}   -input_type  fasta   -dbtype  nucl  -out   {params.GnmBlasdbID}"




rule blast:
	input:
		GnmSimulReadsF="/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.simulateReads.fasta",
		GnmBlasdbF="/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.blastdb.ndb",
	output:
		ReadsBlastOF="/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.blastOut.txt",
	params:
		GnmBlasdbID="/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.blastdb",
	shell:
		"blastn -db    {params.GnmBlasdbID}   -query  {input.GnmSimulReadsF}  -outfmt 6 -num_threads 4   1>{output.ReadsBlastOF}"


rule select_repeatRegion:
	input:
		ReadsBlastOF="/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.blastOut.txt",

	output:
		"/shared/liuhj/HP/process/Genome_repeatRegion/{samp}/{samp}.RepeatRegion.txt"
	params:
		refID="{samp}",
		scriptP="/shared/liuhj/HP/scripts/HP-microevolution/scripts/genome_repeatRegion",
	
	shell:
		"python {params.scriptP}/Select-RepeatRegion-from-readsBlastF.py  -i {input[0]} -o  {output}  -r {params.refID} "






