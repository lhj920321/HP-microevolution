#服务器
necat_path=/home/liuhj/Documents/software/NECAT/Linux-amd64/bin   #docker
ONT_data_path=/home/data/ONT/third   #docker

sample=$1
sample_ONT_FILE=$2
NECAT_out_path=$3
samp_ONT_config_path=$4

echo $sample
echo $sample_ONT_FILE
echo $NECAT_out_path
echo $samp_ONT_config_path
mkdir $samp_ONT_config_path

samp_ONT_READ_LIST=$NECAT_out_path/${sample}.read_list.txt
samp_config_file_1=$samp_ONT_config_path/sample.${sample}.config_file_1.txt
samp_config_file=$samp_ONT_config_path/sample.${sample}.config_file.txt

echo $samp_config_file_1
echo $samp_config_file

for (( step = 0; step  <= 4; step ++ )); do

#step == 0  prepare data fastq_merge
#step == 1  config
#step == 2  correct
#step == 3  assembly
#step == 4  bridge


GENOME_SIZE=1600000
THREADS=80
MIN_READ_LENGTH=5000
CNS_OUTPUT_COVERAGE=40
NUM_ITER=1


OVLP_FAST_OPTIONS_n=500
OVLP_FAST_OPTIONS_z=20
OVLP_FAST_OPTIONS_b=2000
OVLP_FAST_OPTIONS_e=0.5
OVLP_FAST_OPTIONS_j=0
OVLP_FAST_OPTIONS_u=1
OVLP_FAST_OPTIONS_a=1000


OVLP_SENSITIVE_OPTIONS_n=500
OVLP_SENSITIVE_OPTIONS_z=10
OVLP_SENSITIVE_OPTIONS_e=0.5
OVLP_SENSITIVE_OPTIONS_j=0
OVLP_SENSITIVE_OPTIONS_u=1
OVLP_SENSITIVE_OPTIONS_a=1000

CNS_FAST_OPTIONS_a=2000
CNS_FAST_OPTIONS_x=4
CNS_FAST_OPTIONS_y=12
CNS_FAST_OPTIONS_l=1000
CNS_FAST_OPTIONS_e=0.5
CNS_FAST_OPTIONS_p=0.8
CNS_FAST_OPTIONS_u=0

CNS_SENSITIVE_OPTIONS_a=2000
CNS_SENSITIVE_OPTIONS_x=4
CNS_SENSITIVE_OPTIONS_y=12
CNS_SENSITIVE_OPTIONS_l=1000
CNS_SENSITIVE_OPTIONS_e=0.5
CNS_SENSITIVE_OPTIONS_p=0.8
CNS_SENSITIVE_OPTIONS_u=0

TRIM_OVLP_OPTIONS_n=100
TRIM_OVLP_OPTIONS_z=10
TRIM_OVLP_OPTIONS_b=2000
TRIM_OVLP_OPTIONS_e=0.5
TRIM_OVLP_OPTIONS_j=1
TRIM_OVLP_OPTIONS_u=1
TRIM_OVLP_OPTIONS_a=400

ASM_OVLP_OPTIONS_n=100
ASM_OVLP_OPTIONS_z=10
ASM_OVLP_OPTIONS_b=2000
ASM_OVLP_OPTIONS_e=0.5
ASM_OVLP_OPTIONS_j=1
ASM_OVLP_OPTIONS_u=0
ASM_OVLP_OPTIONS_a=400

##prepare data
if [[ $step == 0 ]]; then
    echo $sample_ONT_FILE >$samp_ONT_READ_LIST
fi



if [[ $step == 1 ]]; then
	## config
	#cd $samp_ONT_config_path
    $necat_path/necat.pl  config   $samp_config_file_1

    ##set
    rm $samp_config_file   ##delete raw config file
    cat $samp_config_file_1  | while read line; do
    	#echo $line
    	if [[ $line =~ 'PROJECT' ]]; then
    		out='PROJECT='${sample}

        elif [[ $line =~ 'ONT_READ_LIST' ]]; then
            out='ONT_READ_LIST'=${samp_ONT_READ_LIST}

    	elif [[ $line =~ 'GENOME_SIZE' ]]; then
    		out='GENOME_SIZE='${GENOME_SIZE}

    	elif [[ $line =~ 'THREADS' ]]; then
    		out='THREADS='${THREADS}

    	elif [[ $line =~  'MIN_READ_LENGTH' ]]; then
    		out='MIN_READ_LENGTH='${MIN_READ_LENGTH}
        
        elif [[ $line =~ 'OVLP_FAST_OPTIONS' ]]; then
            out='OVLP_FAST_OPTIONS="-n '${OVLP_FAST_OPTIONS_n}' -z '${OVLP_FAST_OPTIONS_z}' -b '${OVLP_FAST_OPTIONS_b}' \
            -e '${OVLP_FAST_OPTIONS_e}' -j '${OVLP_FAST_OPTIONS_j}' -u '${OVLP_FAST_OPTIONS_u}' -a '${OVLP_FAST_OPTIONS_a}'"'

        elif [[ $line =~ 'OVLP_SENSITIVE_OPTIONS' ]]; then
            out='OVLP_SENSITIVE_OPTIONS="-n '${OVLP_SENSITIVE_OPTIONS_n}' -z '${OVLP_SENSITIVE_OPTIONS_z}' -e '${OVLP_SENSITIVE_OPTIONS_e}' -j '${OVLP_SENSITIVE_OPTIONS_j}' -u '${OVLP_SENSITIVE_OPTIONS_u}' -a '${OVLP_SENSITIVE_OPTIONS_a}'"'

        elif [[ $line =~ 'CNS_FAST_OPTIONS' ]]; then
            out='CNS_FAST_OPTIONS="-a '${CNS_FAST_OPTIONS_a}' -x '${CNS_FAST_OPTIONS_x}' -y '${CNS_FAST_OPTIONS_y}' -l '${CNS_FAST_OPTIONS_l}' -e '${CNS_FAST_OPTIONS_e}' -p '${CNS_FAST_OPTIONS_p}' -u '${CNS_FAST_OPTIONS_u}'"'

        elif [[ $line =~ 'CNS_SENSITIVE_OPTIONS' ]]; then
            out='CNS_SENSITIVE_OPTIONS="-a '${CNS_SENSITIVE_OPTIONS_a}' -x '${CNS_SENSITIVE_OPTIONS_x}' -y '${CNS_SENSITIVE_OPTIONS_y}' -l '${CNS_SENSITIVE_OPTIONS_l}' -e '${CNS_SENSITIVE_OPTIONS_e}' -p '${CNS_SENSITIVE_OPTIONS_p}' -u '${CNS_SENSITIVE_OPTIONS_u}'"'

        elif [[ $line =~ 'TRIM_OVLP_OPTIONS' ]]; then
            out='TRIM_OVLP_OPTIONS="-n '${TRIM_OVLP_OPTIONS_n}' -z '${TRIM_OVLP_OPTIONS_z}' -b '${TRIM_OVLP_OPTIONS_b}' -e '${TRIM_OVLP_OPTIONS_e}' -j '${TRIM_OVLP_OPTIONS_j}' -u '${TRIM_OVLP_OPTIONS_u}' -a '${TRIM_OVLP_OPTIONS_a}'"'

        elif [[ $line =~ 'ASM_OVLP_OPTIONS' ]]; then
            out='ASM_OVLP_OPTIONS="-n '${ASM_OVLP_OPTIONS_n}' -z '${ASM_OVLP_OPTIONS_z}' -b '${ASM_OVLP_OPTIONS_b}' -e '${ASM_OVLP_OPTIONS_e}' -j '${ASM_OVLP_OPTIONS_j}' -u '${ASM_OVLP_OPTIONS_u}' -a '${ASM_OVLP_OPTIONS_a}'"'

        elif [[ $line =~ 'NUM_ITER=2' ]]; then
            out='NUM_ITER='${NUM_ITER}

        elif [[ $line =~ 'CNS_OUTPUT_COVERAGE' ]]; then
            out='CNS_OUTPUT_COVERAGE='${CNS_OUTPUT_COVERAGE}

    	else
    		out=$line
    	fi

    	echo $out >>$samp_config_file
    done
    rm $samp_config_file_1 
fi



if [[ $step == 2 ]]; then
	## correct
	cd $NECAT_out_path
    $necat_path/necat.pl  correct   $samp_config_file
fi


if [[ $step == 3 ]]; then
	## assemble
	cd $NECAT_out_path
    $necat_path/necat.pl  assemble  $samp_config_file
fi

if [[ $step == 4 ]]; then
    ## bridge
    cd $NECAT_out_path
    $necat_path/necat.pl  bridge  $samp_config_file
fi


 
done


#  cp    /home/liuhj/liuhj_2_44/nanopore_HP_1/process/NECAT/pipeline_NECAT_genome/{samp}.bridged_contigs.fasta




