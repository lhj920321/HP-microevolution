#服务器

#1
#BGIseq fastqc   
function fastQC_run(){
   fastqc  -o  $fastqc_output  -f  fastq    $BGIseq_data_R1  $BGIseq_data_R2
}


function pilon_run(){
  echo $all_ONT_NECAT_genome
  echo $pilon_output/${sample}.pilon.0.fasta
  echo $Trim_out_R1_Paired
  echo $Trim_out_R2_Paired
  cp $all_ONT_NECAT_genome   $pilon_output/${sample}.pilon.0.fasta

  for((integer = 1; integer <= ${pilon_times}; integer++))

  do
      echo 
      echo
      echo "item=""$integer"
      echo 
      echo

      ref_index=$(($integer - 1))
      echo 
      echo
      echo $pilon_output/${sample}.${ref_index}.fasta
      echo 
      echo
      bwa index  $pilon_output/${sample}.pilon.${ref_index}.fasta
      echo "hshshshsh"                         
      bwa mem -t ${samtools_threads} $pilon_output/${sample}.pilon.${ref_index}.fasta  \
      $Trim_out_R1_Paired  $Trim_out_R2_Paired | samtools view -bS - > $pilon_output/${sample}.${ref_index}.bam
      samtools sort  --threads  ${samtools_threads}  -o $pilon_output/${sample}.${ref_index}.sort.bam  $pilon_output/${sample}.${ref_index}.bam
      samtools index $pilon_output/${sample}.${ref_index}.sort.bam
 

      java   -Xmx60g  -jar  /home/amax/anaconda3/envs/biotools_lhj/share/pilon-1.24-0/pilon.jar --genome  $pilon_output/${sample}.pilon.${ref_index}.fasta \
      --bam $pilon_output/${sample}.${ref_index}.sort.bam \
      --output ${sample}.pilon.${integer} \
      --outdir $pilon_output --threads ${samtools_threads}   # 2>> $pilon_output/${sample}.pilon.log   

#      java -jar  /home/amax/anaconda3/envs/biotools_lhj/share/pilon-1.24-0/pilon-1.24.jar 
      #--threads ${samtools_threads}
      #rm $pilon_output/${sample}.${ref_index}.sort.bam

      echo ${sample}.pilon.${integer}
  done
}


##################################################


#set input 参数
sample=$1
run_step=$2
samtools_threads=30

#java_1_8_path=/home/liuhj/liuhj_2_44/software/jre1.8.0_191/bin
#pilon_path=/home/liuhj/liuhj_2_44/nanopore_HP_1/software




if [[ $run_step == fastQC ]]; then
##fastqc
BGIseq_data_R1=$3
BGIseq_data_R2=$4
##output
fastqc_output=$5


##shell
	fastQC_run
fi



if [[ $run_step == pilon ]]; then

#input
  all_ONT_NECAT_genome=$3
  Trim_out_R1_Paired=$4  
  Trim_out_R2_Paired=$5  

##output
  pilon_output=$6
  pilon_times=$7

#input
  #all_ONT_NECAT_genome=/public1/home/liuhj2/Nanopore_HP_2/process/NECAT/${sample}.necat.fasta
  #Trim_out_R1_Paired=/public1/home/liuhj2/Nanopore_HP_2/process/Trimmomatic/${sample}/${sample}.R1.paired.fastq  
 # Trim_out_R2_Paired=/public1/home/liuhj2/Nanopore_HP_2/process/Trimmomatic/${sample}/${sample}.R2.paired.fastq  

##output
  #pilon_output=/public1/home/liuhj2/Nanopore_HP_2/process/pilon/${sample}


  rm -r $pilon_output
  mkdir $pilon_output
##shell
  pilon_run


fi









