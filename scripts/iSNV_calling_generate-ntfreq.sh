#!/bin/bash

DEP_THRES=50    ##最低深度
MIN_DEP_THRES=5   ##支持alt的最小reads数
VALID_SIZE_THRES=100
FREQ_THRES=0.02   ###alt的频率最小值
STRANDED_RATIO_THRES=0.1     ## reads链偏性检测
# #download sickle and hammer/spades, install the executables in your $PATH
#sickleBinPath="~/bin/sickle-master"
#spadesBinPath="~/bin/SPAdes-3.7.1-Linux/"/public1/home/liuhj/project_2
fastqPath="./sample_reads"
cleanReadPath=/public1/home/liuhj2/process/Trimmomatic/${sample}
noThreadStep1=11
samtoolsThread=4
bowtieThread=24


customStep=$1   ##3,4,5
#ref_sample=$2

#refFasta=/public1/home/liuhj2/process/pilon/${ref_samp}/${ref_samp}.piloned.fasta     #"./sample_ref/Ebola_genome1.fa"
bowtie2indexPath=/public1/home/liuhj2/process/bwa/index/BGI_bwa2All_ONT_NECAT #"./sample_ref/Ebola_genome"
mpileupPath=/shared/liuhj/coffee/tonglv/Liangci/ntfreq/Same_sample/0147_0148             #"./mpileup_and_ntfreq"

noThreadStep3=6
tablePath=/shared/liuhj/tonglv/process/raw_ntfreq/table_0147_0148
mkdir $tablePath

excludeRegionF=/shared/liuhj/tonglv/process/ntfreq/exclude_region.txt


refFasta="/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna"     #"./sample_ref/Ebola_genome1.fa"



#########################################################################3
echo 
echo "********************************************************"
echo "*       Intrahost SNV calling for amplicon-seq         *"
echo "********************************************************"
echo "nim, 201709 v2"
echo


if [[ $customStep =~ "ntfreq" ]];then
bash step3_mpileup2ntfreq.sh $noThreadStep3 $mpileupPath



