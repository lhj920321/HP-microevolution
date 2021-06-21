#!/bin/bash

## 
AAfileP=$1
examp_ctlF=$2
out_log=$3

echo "control file "



for AAfile in $AAfileP/*.aa.fna; do
	CodonFile=`echo ${AAfile%.aa*}.codon | sed "s/\//--/g"`
	#TreeFile=`echo ${AAfile%.aa*}.treefile | sed "s/\//__/g"`
	codemlOutF=`echo ${AAfile%.aa*}.codeml | sed "s/\//--/g"`

	outF=${AAfile%.aa*}.ctl


	cp  $examp_ctlF  $outF
	sed  -i   "s/CodonFile/$CodonFile/" $outF
	#sed  -i   "s/TreeFile/$TreeFile/"   $outF
	sed  -i   "s/codemlOutF/$codemlOutF/" $outF




	sed  -i  "s/--/\//g"   $outF

	echo $CodonFile
	echo $outF




done



echo "" >$out_log




