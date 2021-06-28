###################################
#######20210623
#######phasing in iSNV peak  based on the NGS and CCS reliable posi and CCS alignment
###################################
#minFreq="0.86"
#maxFreq="0.96"

#需要根据每个样本不同的峰区域 进行调整
CCS_Samps=["P6-E-d"]
Ref_samp="P6-E-t"


rule all:
	input:
	##read_LongReadbam_file
	#0 filt SNV posi
		expand("/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-Filted_reliablePosi.txt",samp=CCS_Samps,Ref_samp=Ref_samp),
		#expand("/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-Filted_reliablePosi.txt",samp=CCS_Samps,Ref_samp=Ref_samp),
	#1 read bam -reads haplotype
		expand("/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-SNVtype_of_Filted_Cites.txt",samp=CCS_Samps,Ref_samp=Ref_samp),
	#2 phasing 
		expand("/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/step2_phasing.log",samp=CCS_Samps,Ref_samp=Ref_samp),

#0
rule FiltNGS_CCS_Posi:
	input:
		CCS_VcfF = "/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/table/allVcf/{samp}_2_{Ref_samp}.minFreq0.02.iSNVpy.vcf",
		NGS_VcfF = "/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/allVcf/{samp}_2_{Ref_samp}.minFreq0.02.iSNVpy.vcf",
	output:
		out_reliableCite_CCS = "/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-Filted_reliablePosi.txt",
		#out_0_line = "/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.log.txt",
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/phasing-Gnm/CCS",
		minFreq="0.86",
		maxFreq="0.96",
		minFreqReliaPosi="0.02"
	shell:
		"python {params.scriptPath}/Step0__read_vcf_FiltPosi.py  -N {input.NGS_VcfF} \
		-P {input.CCS_VcfF} -r {params.minFreqReliaPosi} -m   {params.minFreq} -M   {params.maxFreq} -o {output.out_reliableCite_CCS} "



#1
rule read_LongReadbam_file:
	input:
		configF="/shared/liuhj/HP/process/phasing_CCS/config.txt",
		reliablePosiVcf="/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-Filted_reliablePosi.txt",
		BamF="/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/ngmlr/{samp}/{samp}_2_{Ref_samp}.sort.bam",
		refSampGnm="/shared/liuhj/HP/process/assembly/all_genomes/{Ref_samp}.fasta",
	output:	
		SNV_type="/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-SNVtype_of_Filted_Cites.txt",
		SNV_Type_Stat="/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-SNVtype_of_Filted_Cites.Stat.txt",
		#"/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-bam_SnvType.txt",
		#"/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-Filted_Cites.txt",
		#"/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-Filted_reliablePosi.txt",		
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/phasing-Gnm/CCS",
		outP="/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96"
	shell:
		"python {params.scriptPath}/Step1__read_ONTbam_f.py  -c {input.configF}   -b {input.BamF} -v {input.reliablePosiVcf} -r {input.refSampGnm} -o {params.outP}"  


#2
rule phasing:
	input:
		reliablePosiVcf="/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-Filted_reliablePosi.txt",
		SNV_Type_Stat="/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/1.out-SNVtype_of_Filted_Cites.Stat.txt",
	output:
		phasingLog="/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96/step2_phasing.log",
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/phasing-Gnm/CCS-1",
		outP="/shared/liuhj/HP/process/phasing_CCS-1/phasing_{samp}_2_{Ref_samp}--0.86-0.96",
		startLociNum="10",
		minLociNum_forOut="20",
		minCov="0.7",

	shell:
		"python {params.scriptPath}/step3_phasing.py  -f {input.reliablePosiVcf} -F {input.SNV_Type_Stat} -o {params.outP}  -s {params.startLociNum} -m {params.minLociNum_forOut} -c {params.minCov} -l {output.phasingLog}"





