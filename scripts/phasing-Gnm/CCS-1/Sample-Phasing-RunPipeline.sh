#


SampConfigF="/shared/liuhj/HP/scripts/HP-microevolution/scripts/phasing-Gnm/CCS-1/Sample-config.txt"
Snakemake_templateF="/shared/liuhj/HP/scripts/HP-microevolution/scripts/phasing-Gnm/CCS-1/phasing.snakemake.phasing.py"
Samp_snakemakeP="/shared/liuhj/HP/scripts/HP-microevolution/scripts/phasing-Gnm/CCS-1/Samp_phasing_RunScripts"
mkdir $Samp_snakemakeP

cat $SampConfigF | while read line ; do
	echo $line
	RUN_Sample=`echo $line | awk '{print $1}'`
	RefSample_MapTo=`echo $line | awk  '{print $2}'`
	MIN_Freq_Of_feng=`echo $line | awk '{print $3}'`
	MAX_Freq_Of_feng=`echo $line | awk '{print $4}'`

	echo $RUN_Sample
SampScriptF=$Samp_snakemakeP/${RUN_Sample}--${MIN_Freq_Of_feng}-${MAX_Freq_Of_feng}.phasing.snakemake.phasing.py
##cp 
#if [ ! -f  $SampScriptF ]; then 
echo "Sample snakemake script not exist!" 
cp $Snakemake_templateF  $SampScriptF
echo $SampConfigF
echo $SampScriptF

echo ${RUN_Sample}
sed -i "s/RUN_Sample/${RUN_Sample}/"  $SampScriptF
sed -i "s/RefSample_MapTo/${RefSample_MapTo}/"  $SampScriptF
sed -i "s/MIN_Freq_Of_feng/${MIN_Freq_Of_feng}/"  $SampScriptF
sed -i "s/MAX_Freq_Of_feng/${MAX_Freq_Of_feng}/"  $SampScriptF




snakemake -s $SampScriptF -p  -j 1
#snakemake -s $SampScriptF -np


done







