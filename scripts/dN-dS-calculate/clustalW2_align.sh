#!/bin/bash

## align the AA sequences
AAfileP=$1
out_log=$2
echo "clustalw2 align"
for AAfile in $AAfileP/*.aa.fna; do
	echo $AAfile
	outalnF=${AAfile%.fna*}.aln
	outdndF=${AAfile%.fna*}.dnd

	echo $outalnF
	echo -e "1\n${AAfile}\n2\n1\n${outalnF}\n${outdndF}\nX\n\nX\n"  | clustalw2

	rm $outdndF

done
echo "" >$out_log






