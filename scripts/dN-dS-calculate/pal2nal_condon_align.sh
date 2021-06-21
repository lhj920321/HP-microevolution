#!/bin/bash

## pal2nal codon file
AAfileP=$1
out_log=$2
echo "pal2nal  codon"
for AAfile in $AAfileP/*.aa.fna; do
	echo $AAfile
	outalnF=${AAfile%.fna*}.aln
	NucF=${AAfile%.aa*}.nuc.fna
	outF=${AAfile%.aa*}.codon

	echo $outalnF
	echo $NucF
	echo $outF

	pal2nal.pl  $outalnF  $NucF -output paml   -nogap  >$outF

done




echo "" >$out_log










