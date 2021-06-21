#虚拟机

person=$1


if [[ $person == 'P1' ]]; then
	ref_samp='P1-E-j'
	samp_list="P1-E-x P1-E-d P1-E-j P1-E-t"
	#samp_list="P1-E-x"
fi

if [[ $person == 'P2' ]]; then
	ref_samp="P2-E-t"
	samp_list="P2-E-x P2-E-d P2-E-j P2-E-t"
fi

if [[ $person == 'P3' ]]; then
	ref_samp="P3-E-j"
	samp_list="P3-E-x P3-E-j P3-E-t"
fi

if [[ $person == 'P4' ]]; then
	ref_samp="P4-E-j"
	samp_list="P4-E-x P4-E-d P4-E-j"
fi

if [[ $person == 'P5' ]]; then
	ref_samp="P5-E-d"
	samp_list="P5-E-d P5-E-j"
fi

if [[ $person == 'P6' ]]; then
	ref_samp="P6-E-t"
	samp_list="P6-E-d P6-E-t"
fi

if [[ $person == 'P7' ]]; then
	ref_samp="P7-E-j"
	samp_list="P7-E-d P7-E-j P7-E-t"
fi


if [[ $person == 'P8' ]]; then
	ref_samp="P8-E-x"
	samp_list="P8-E-x P8-E-d P8-E-t"
fi

if [[ $person == 'P9' ]]; then
	ref_samp="P9-E-j"
	samp_list="P9-E-d P9-E-j P9-E-t"
fi

if [[ $person == 'P10' ]]; then
	ref_samp="P10-E-j"
	samp_list="P10-E-x P10-E-j P10-E-t"
fi


if [[ $person == 'P11' ]]; then
	ref_samp="P11-E-x"
	#samp_list="P11-E-x P11-E-j P11-E-t"
	samp_list="P11-E-x P11-E-t"
fi


if [[ $person == 'P12' ]]; then
	ref_samp="P12-E-t"
	samp_list="P12-E-x P12-E-t"
fi

if [[ $person == 'P13' ]]; then
	ref_samp="P13-E-j"
	samp_list="P13-E-x P13-E-d P13-E-j P13-E-t"
fi

if [[ $person == 'P14' ]]; then
	ref_samp="P14-E-j"
	samp_list="P14-E-x P14-E-d P14-E-j"
fi

if [[ $person == 'P15' ]]; then
	ref_samp="P15-C_bS-j"
	samp_list="P15-C_bS-d P15-C_bS-j P15-C_bS-t"
fi

if [[ $person == 'P16' ]]; then
	ref_samp="P16-C_bS-j"
	samp_list="P16-C_bS-x P16-C_bS-d P16-C_bS-j"
fi

if [[ $person == 'P17' ]]; then
	ref_samp="P17-C_bR-j"
	samp_list="P17-C_bR-x P17-C_bR-d P17-C_bR-j"
fi

if [[ $person == 'P18' ]]; then
	ref_samp="P18-C_bR-x"
	samp_list="P18-C_bR-x P18-C_bR-j"
fi

if [[ $person == 'P19' ]]; then
	ref_samp="P19-C_bR-j"
	samp_list="P19-C_bR-d P19-C_bR-j P19-C_bR-t"
fi

if [[ $person == 'P20' ]]; then
	ref_samp="P20-C_bR-t"  ##changed
	samp_list="P20-C_bR-x P20-C_bR-d P20-C_bR-j P20-C_bR-t"
fi

if [[ $person == 'P21' ]]; then
	ref_samp="P21-C_sL-j"
	samp_list="P21-C_sL-x P21-C_sL-d P21-C_sL-j P21-C_sL-t"
fi

if [[ $person == 'P22' ]]; then
	ref_samp="P22-C_sL-d"  ##changed
	samp_list="P22-C_sL-d P22-C_sL-t"
fi

if [[ $person == 'P23' ]]; then
	ref_samp="P23-C_sL-j"
	samp_list="P23-C_sL-x P23-C_sL-j"
fi

if [[ $person == 'P24' ]]; then
	ref_samp="P24-C_sC-t"  ##changed
	samp_list="P24-C_sC-x P24-C_sC-j P24-C_sC-t"
fi

if [[ $person == 'P25' ]]; then
	ref_samp="P25-C_sC-t"
	samp_list="P25-C_sC-j P25-C_sC-t"
fi


echo $samp_list
echo $ref_samp


for samp in $samp_list; do

#####snpEff
snpEff_P=/shared/liuhj/software/snpEff
snpEff=${snpEff_P}/snpEff.jar
snpEff_config_F=${snpEff_P}/snpEff.config
prokka_out_dir=/shared/liuhj/HP/process/assembly/prokka/${ref_samp}_prokka

Anno_vcfP=/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/${person}/vcf-Anno
mkdir $Anno_vcfP


vcf_P=/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/${person}/vcf
vcf_F=$vcf_P/${samp}_2_${ref_samp}.minFreq0.02.iSNVpy.vcf

#snpEff config
echo $ref_samp
ref_ID=HP_${ref_samp}
hang_1=`echo -e "# HP genome,version  HP_${ref_samp}"`
hang_2=`echo -e "HP_${ref_samp}.genome : HP_${ref_samp}"`

echo $hang_1
echo $hang_2


sed -i "54 a ${hang_1}" $snpEff_config_F
sed -i "55 a ${hang_2}" $snpEff_config_F
databases_gff_P=${snpEff_P}/data/HP_${ref_samp}
databases_fna_P=${snpEff_P}/data/genome

mkdir  $databases_gff_P

gff_F=$prokka_out_dir/${ref_samp}.gff
fna_F=$prokka_out_dir/${ref_samp}.fna

databases_gff_F=$databases_gff_P/${ref_samp}.gff
databases_fna_F=$databases_fna_P/${ref_samp}.fna
cp  $gff_F  $databases_gff_F
cp  $fna_F  $databases_fna_F
mv $databases_gff_F  $databases_gff_P/genes.gff
mv $databases_fna_F  $databases_fna_P/${ref_ID}.fna
cd  $snpEff_P
java -jar  $snpEff  build  -gff3  -v ${ref_ID}

touch $Anno_vcfP/${samp}_2_${ref_samp}.minFreq0.02.snpEffAnno.vcf
java -jar  $snpEff   $ref_ID   $vcf_F   >$Anno_vcfP/${samp}_2_${ref_samp}.minFreq0.02.snpEffAnno.vcf
mv $snpEff_P/snpEff_summary.html  $Anno_vcfP/${samp}.minFreq0.02.snpEff.summary.html



done




