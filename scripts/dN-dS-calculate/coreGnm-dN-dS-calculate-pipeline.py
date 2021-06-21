###################################
###################################
#######CCS data map to respe genome
###################################


personLst = ["P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P12","P15","P16","P17","P18","P19","P20","P21","P22","P23","P24","P25"]

#"P13","P14","P11","P8"
personLst = ["P13","P14"]

rule all:
	input:
	## nuc and aa file (diff gene)
		expand("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}",person=personLst),
	##clustalW
		#expand("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.clustalW.log",person=personLst),
	##pal2nal codon align
		#expand("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.pal2nal.log",person=personLst),
	## control file 
		#expand("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.ctl.log",person=personLst),
	##codeml run (dN dS)
		#expand("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.codeml.log",person=personLst),
	## control file  -null
		#expand("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.ctl.null.log",person=personLst),
	##codeml run (dN dS) -null
		#expand("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.codeml.null.log",person=personLst),

	##codeml secect dN dS
		expand("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select/{person}.codeml.dNdS.txt",person=personLst),



'''
rule Nuc_AA_FnaF_of_diverGene:
	input:	
		roaryP = "/shared/liuhj/HP/process/samps_Tree-20210527",
	output:
		directory("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}"),
	params:
		prokkaP = "/shared/liuhj/HP/process/assembly/prokka",
		scriptP = "/shared/liuhj/HP/scripts/HP-microevolution/scripts/dN-dS-calculate",
		personID = "{person}",
	shell:
		"python {params.scriptP}/Core_Genome_nuc_aa-seq-from-roaryOutput.py  -i {params.personID} -r {input.roaryP}  -p {params.prokkaP}  -o  {output}"  
'''



rule clustalW2_align:
	input:
		"/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}"
	output:
		temp("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.clustalW.log")
	params:
		scriptP = "/shared/liuhj/HP/scripts/HP-microevolution/scripts/dN-dS-calculate",
	shell:
		"bash {params.scriptP}/clustalW2_align.sh  {input}  {output}"






rule pal2nal_codon:
	input:
		"/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}",
		#temp("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.clustalW.log"),
	output:
		"/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.pal2nal.log"
	params:
		scriptP = "/shared/liuhj/HP/scripts/HP-microevolution/scripts/dN-dS-calculate",
	shell:
		"bash {params.scriptP}/pal2nal_condon_align.sh  {input[0]}  {output}"





rule Generate_controlF:
	input:
		PersonP="/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}",
		exampF="/shared/liuhj/HP/process/dN-dS-caculate/paml-control-file-sample.ctl",
		pal2nalLog=temp("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.pal2nal.log")
	output:
		"/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.ctl.log",
	params:
		scriptP = "/shared/liuhj/HP/scripts/HP-microevolution/scripts/dN-dS-calculate",
	shell:
		"bash {params.scriptP}/Generate_sample_controlF.sh  {input.PersonP}  {input.exampF}  {output}"





rule codeml_run:
	input:
		ctlLog=temp("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.ctl.log"),
		PersonP="/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}",
	output:
		temp("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.codeml.log"),

	params:
		scriptP = "/shared/liuhj/HP/scripts/HP-microevolution/scripts/dN-dS-calculate",
		personID="{person}",
	shell:
		"bash {params.scriptP}/codeml-run.sh  {input.PersonP} {output}"



##Null
'''
rule Generate_controlF_null:
	input:
		PersonP="/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}",
		exampF="/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/paml-control-file-sample-null.ctl",
		#pal2nalLog=temp("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.pal2nal.log")
	output:
		"/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.ctl.null.log",
	params:
		scriptP = "/shared/liuhj/HP/scripts/HP-microevolution/scripts/dN-dS-calculate",
	shell:
		"bash {params.scriptP}/Generate_sample_controlF-null.sh  {input.PersonP}  {input.exampF}  {output}"






rule codeml_run_null:
	input:
		ctlLog=temp("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.ctl.null.log"),
		PersonP="/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}",
	output:
		temp("/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.codeml.null.log"),

	params:
		scriptP = "/shared/liuhj/HP/scripts/HP-microevolution/scripts/dN-dS-calculate",
		personID="{person}",
	shell:
		"bash {params.scriptP}/codeml-run-null.sh  {input.PersonP} {output}"
'''




rule codeml_dNdS_select:
	input:
		PersonP="/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}",
		#codemlLog="/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/{person}.codeml.log"
	output:
		"/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select/{person}.codeml.dNdS.txt",

	params:
		scriptP = "/shared/liuhj/HP/scripts/HP-microevolution/scripts/dN-dS-calculate",

	shell:
		"python {params.scriptP}/colect_dN_dS.py  -c {input.PersonP}  -o  {output} "







