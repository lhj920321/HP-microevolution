#!/bin/bash

## codeml run for dN and dS
AAfileP=$1
out_log=$2


echo "codeml run for dN and dS "

for AAfile in $AAfileP/*.aa.fna; do

	ctlF=${AAfile%.aa*}.null.ctl

	echo "control file : "$ctlF
	codeml $ctlF



done



echo "" >$out_log










