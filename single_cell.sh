#!/bin/bash

# qsub -W group_list=cge -A cge -l nodes=1:ppn=40,mem=50g,walltime=40:00:00 sc.queue.sh

################################################################################################################
############################################## help page #######################################################
Help()
{
   # Display Help
   echo -e "backup your files first !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
   echo "Add your parameters here."
   echo "options:"
   echo "-i|--input   Necessary. Path to the raw undemultiplexed fastq/fastq.gz files directory. Sample names have to be named as by 'sample_R1*' or sample_R2*."
   echo "-r|--reference  Necessary. Path to the reference directory. Sample names have to be seperated by '_'."
   echo "-o|--output    Necessary. Path to the output directory"
   echo "--barcode1 Necessary. The barcode1.fasta file path. The barcode1 is the one closest to the reads end."
   echo "--barcode2 Necessary. The barcode2.fasta file path. The barcode2 is the one 2nd closest to the reads end."
   echo "--barcode3 Necessary. The barcode3.fasta file path. The barcode3 is the one 3rd closest to the reads end."
   echo "--primer  Necessary. The primer fasta file path."
   echo "--minimal_reads_number   Optional. The minimal number of reads R1+R2, a value or "". "
   echo "--maximal_reads_number   Optional. The maximal number of reads R1+R2, a value or "". "
   echo "-t|--threads   Necessary. Number of threads to use"   
   echo "--trimmomatic_quality  Optional. The quality score for trimmomatic windows."
   echo "--memory   If use "--memory", you will specify the memory to be used (GB)"   
   echo "--cutadapt  Necessary. The executable cutadapt path."
   echo "--trimmomatic  Necessary. The executable trimmomatic path."
   echo "--bwa_mem2  Necessary. The executable bwa-mem2 path."
   echo "--samtools  Necessary. The executable samtools path."
   echo "--bioawk  Necessary. The executable bioawk path."

}

################################################################################################################
############################################## set parameters (for developer) ##################################

# initialize parameters
parameter_memory=20 #in GB
parameter_t=40
parameter_trimmomatic_quality=15

##
VALID_ARGS=$(getopt -o hi:r:o:t: --long input:reference:,output:,threads:,barcode1:,barcode2:,barcode3:,primer:,memory:,\
trimmomatic_quality:,gzip:,cutadapt:,trimmomatic:,bwa_mem2:,samtools:,bioawk:,minimal_reads_number:,maximal_reads_number: -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -h|--help)
     Help
     exit
     ;; 
    -i|--input)
    parameter_i=$2
    shift 2 
    ;;
    -r|--reference)
    parameter_r=$2
    shift 2
    ;;
    -o|--output)
    parameter_o=$2
    shift 2
    ;;
    --barcode1)
    parameter_b1=$2
    shift 2
    ;;
    --barcode2)
    parameter_b2=$2
    shift 2
    ;;
    --barcode3)
    parameter_b3=$2
    shift 2
    ;;
    --primer)
    parameter_primer=$2
    shift 2
    ;;
    -t|--threads)
    parameter_t=$2
    shift 2
    ;;
    --memory)
    parameter_memory=$2
    shift 2
    ;;
    --trimmomatic_quality)
    parameter_trimmomatic_quality=$2
    shift 2
    ;;
    --cutadapt)
    cutadapt_path=$2
    shift 2
    ;;
    --trimmomatic)
    trimmomatic_path=$2
    shift 2
    ;;
    --bwa_mem2)
    bwa_mem2_path=$2
    shift 2
    ;;
    --samtools)
    samtools_path=$2
    shift 2
    ;;
    --bioawk)
    bioawk_path=$2
    shift 2
    ;;
    --minimal_reads_number)
    minimal_reads_number=$2
    shift 2
    ;;
    --maximal_reads_number)
    maximal_reads_number=$2
    shift 2
    ;;
    --) 
    shift;
    break;;
  esac
done

############################################## check parameters 
if [[ ! -z $parameter_i ]]; then echo -e "\nYour input fastq is in $parameter_i"; else echo -e "\nYour input fastq is missing !!!"; echo "use -h for help"; exit;fi
if [[ ! -d $parameter_i ]]; then echo -e "\n$parameter_i does not exist !!!"; exit;fi
if [[ ! -z $parameter_r ]]; then echo -e "\nYour reference for mapping is $parameter_r"; else echo -e "\nYour reference for mapping is missing !!!"; echo "use -h for help"; exit;fi
if [[ ! -d $parameter_r ]]; then echo -e "\n$parameter_r directory does not exit !!!, or you should specify a directory path rather than a file path"; exit; fi
if [ ! -z $parameter_o ]; then echo -e "\nYour output directory is $parameter_o"; else echo -e "\nYour output directory is missing !!!"; echo "use -h for help"; exit;fi
if [[ ! -d $parameter_o ]]; then echo -e "\n$parameter_o does not exist !!!"; exit;fi
if [[ ! -z $parameter_b1 && -f $parameter_b1 ]]; then echo -e "\nYour barcode file 1 is $parameter_b1"; else echo -e "\nYour barcode1 file is missing !!!"; echo "use -h for help"; exit; fi
if [[ ! -z $parameter_b2 && -f $parameter_b2 ]]; then echo -e "\nYour barcode file 2 is $parameter_b2"; else echo -e "\nYour barcode2 file is missing !!!"; echo "use -h for help"; exit; fi
if [[ ! -z $parameter_b3 && -f $parameter_b3 ]]; then echo -e "\nYour barcode file 3 is $parameter_b3"; else echo -e "\nYour barcode3 file is missing !!!"; echo "use -h for help"; exit; fi
if [[ ! -z $parameter_primer && -f $parameter_primer ]]; then echo -e "\nYour primer file is $parameter_primer"; else echo -e "\nYour primer fasta file is missing"; exit; fi

if [[ -z $minimal_reads_number && -z $maximal_reads_number ]]; then
    echo -e "\nYou have to specify one or both of --minimal_reads_number and --maximal_reads_number"
    exit
elif [[ -n $minimal_reads_number && -n $maximal_reads_number ]]; then
    echo -e "\nBoth minimal and maximal reads are specified:"
    echo -e "Minimal reads: $minimal_reads_number"
    echo -e "Maximal reads: $maximal_reads_number"
elif [[ -n $minimal_reads_number ]]; then
    echo -e "\nYour minimal reads is $minimal_reads_number"
elif [[ -n $maximal_reads_number ]]; then
    echo -e "\nYour maximal reads is $maximal_reads_number"
fi

if [[ ! -z $cutadapt_path ]]; then echo -e "\nYour cutadapt is in $cutadapt_path"; else echo -e "\nYour cutadapt is missing !!!"; echo "use -h for help"; exit;fi
if [[ ! -z $trimmomatic_path ]]; then echo -e "\nYour trimmomatic is in $trimmomatic_path"; else echo -e "\nYour trimmomatic is missing !!!"; echo "use -h for help"; exit;fi
if [[ ! -z $samtools_path ]]; then echo -e "\nYour samtools is in $samtools_path"; else echo -e "\nYour samtools is missing !!!"; echo "use -h for help"; exit;fi
if [[ ! -z $bwa_mem2_path ]]; then echo -e "\nYour bwa-mem2 is in $bwa_mem2_path"; else echo -e "\nYour bwa-mem2 is missing !!!"; echo "use -h for help"; exit;fi
if [[ ! -z $bioawk_path ]]; then echo -e "\nYour bioawk is in $bioawk_path"; else echo -e "\nYour bioawk is missing !!!"; echo "use -h for help"; exit;fi
echo "You will use $parameter_t threads"
echo "You will use $parameter_memory GB memory"
echo "You will use trimmomatic quality $parameter_trimmomatic_quality"

################################################################################################################
############################################## file renaming  ##################################################
echo -e "\n Renaming files\n"
ls $parameter_i|while read -r line;do 
new_name1=$(echo $line|sed -e 's/\.R1.fastq/_R1.fastq/' -e 's/\.R1.fq/_R1.fastq/' -e 's/\.1.fastq/_R1.fastq/' -e 's/\.1.fq/_R1.fastq/');\
new_name2=$(echo $new_name1|sed -e 's/\.R2.fastq/_R2.fastq/' -e 's/\.R2.fq/_R2.fastq/' -e 's/\.2.fastq/_R2.fastq/' -e 's/\.2.fq/_R2.fastq/');\
if [[ $line != $new_name2 ]];then 
mv $parameter_i/$line $parameter_i/$new_name2;\
fi;
done 
echo -e "\n Renaming files finished\n"

################################################################################################################
############################################## demultiplexing  ##################################################
echo -e "\n Demultiplexing\n"

out_demtx=$(echo "$parameter_i"|sed 's|$|_demultiplexed|')

# check if demultiplexed dir already exists 
if [[ ! -d $out_demtx ]]; then 
  mkdir -p $out_demtx
  mkdir -p $out_demtx/barcode1
  mkdir -p $out_demtx/barcode2
  mkdir -p $out_demtx/barcode3


  parallel_cutadapt1(){
    fastq1=$1
    fastq2=$(echo $fastq1|sed 's|R1.fastq|R2.fastq|')
    sample=$(basename $fastq1|cut -d_ -f1)
    barcode_name=$(basename $parallel_parameter_b1|cut -d. -f1)
    $parallel_cutadapt_path --action trim -e 0 --cores 1 -g file:$parallel_parameter_b1 -o $parallel_out_demtx/barcode1/${sample}-{name}_R2.fastq -p $parallel_out_demtx/barcode1/${sample}-{name}_R1.fastq $fastq2 $fastq1 >> $parallel_out_demtx/barcode1/${barcode_name}.log 2>&1  
    }
  export parallel_parameter_b1=$parameter_b1
  export parallel_cutadapt_path=$cutadapt_path 
  export parallel_out_demtx=$out_demtx
  export -f parallel_cutadapt1
  parallel -j $parameter_t parallel_cutadapt1 ::: $parameter_i/*R1*
  
  parallel_cutadapt2(){
    fastq1=$1
    fastq2=$(echo $fastq1|sed 's|R1.fastq|R2.fastq|')
    sample=$(basename $fastq1|cut -d_ -f1)
    barcode_name=$(basename $parallel_parameter_b2|cut -d. -f1)
    $parallel_cutadapt_path --action trim -e 0 --cores 1 -g file:$parallel_parameter_b2 -o $parallel_out_demtx/barcode2/${sample}+{name}_R2.fastq -p $parallel_out_demtx/barcode2/${sample}+{name}_R1.fastq $fastq2 $fastq1 >> $parallel_out_demtx/barcode2/${barcode_name}.log 2>&1  
    }
  export parallel_parameter_b2=$parameter_b2
  export parallel_cutadapt_path=$cutadapt_path 
  export parallel_out_demtx=$out_demtx
  export -f parallel_cutadapt2
  ls -lhS $out_demtx/barcode1/*|grep 'R1.fastq'|grep -v 'unknown'|awk '$5 != 0 {print $9}' > $out_demtx/demultiplexing.barcode1.fastq # to skip zero reads files
  parallel -j $parameter_t parallel_cutadapt2 :::: $out_demtx/demultiplexing.barcode1.fastq

  parallel_cutadapt3(){
    fastq1=$1
    fastq2=$(echo $fastq1|sed 's|R1.fastq|R2.fastq|')
    sample=$(basename $fastq1|cut -d_ -f1)
    barcode_name=$(basename $parallel_parameter_b3|cut -d. -f1)
    $parallel_cutadapt_path --action trim -e 0 --cores 1 -g file:$parallel_parameter_b3 -o $parallel_out_demtx/barcode3/${sample}+{name}_R2.fastq -p $parallel_out_demtx/barcode3/${sample}+{name}_R1.fastq $fastq2 $fastq1 >> $parallel_out_demtx/barcode3/${barcode_name}.log 2>&1  
    }
  export parallel_parameter_b3=$parameter_b3
  export parallel_cutadapt_path=$cutadapt_path 
  export parallel_out_demtx=$out_demtx
  export -f parallel_cutadapt3
  ls -lhS $out_demtx/barcode2/*|grep 'R1.fastq'| grep -v 'unknown'|awk '$5 != 0 {print $9}' > $out_demtx/demultiplexing.barcode2.fastq # to skip zero reads files
  parallel -j $parameter_t parallel_cutadapt3 :::: $out_demtx/demultiplexing.barcode2.fastq

  echo -e "\n Demultiplexing finished\n"
fi

################################################################################################################
############################################## remove samples based on reads number  ###########################
echo -e "\n removing samples based on reads number"

if [[ ! -f $out_demtx/fastq.non.zero.reads ]]; then 
    # prune based on sum of R1 + R2
    ls -lhS "$out_demtx/barcode3"|grep -v 'unknown'|awk '$5 != 0 {print $9}'|sed '1d' > $out_demtx/fastq.non.zero # step 1 prune based on file size, since 50% of files are empty by file size 

    process_file() {
    input=$1
    reads_num=$("$parallel_bioawk_path" -cfastx '{if (length($seq)>0) sum=sum+1}END{print sum}' $input) #we count reads when length > 0, sometimes it can have reads name but not sequences so grep -c '^@' will always give same R1 R2 reads number, but our method not necessary R1 always equal R2 number. 
    if [[ "$reads_num" -gt 0 ]]; then 
    output=$(basename "$input")
    echo $output $reads_num
    fi
    }
    export -f process_file # in process_file we keep all reads > 0
    export parallel_bioawk_path="$bioawk_path"
    sed "s|^|$out_demtx/barcode3/|" $out_demtx/fastq.non.zero | parallel -j $parameter_t process_file >> $out_demtx/fastq.non.zero.reads # calculate reads
    awk '{split($1,arr,"_");print arr[1],$0}' $out_demtx/fastq.non.zero.reads > $out_demtx/tmp && mv $out_demtx/tmp $out_demtx/fastq.non.zero.reads
fi

if [[ -n $minimal_reads_number && -n $maximal_reads_number ]]; then 
        awk -v var1="$minimal_reads_number" -v var2="$maximal_reads_number" '{arr[$1] += $3} END {for (key in arr) {if (arr[key]>=var1 && arr[key]<=var2) print key}}' $out_demtx/fastq.non.zero.reads > $out_demtx/filtered.fastq.name # step 2 prune based on reads R1+R2
    elif [[ -n $minimal_reads_number ]]; then 
        awk -v var="$minimal_reads_number" '{arr[$1] += $3} END {for (key in arr) {if (arr[key]>=var) print key}}' $out_demtx/fastq.non.zero.reads > $out_demtx/filtered.fastq.name # step 2 prune based on reads R1+R2
    elif [[ -n $maximal_reads_number ]]; then 
        awk -v var="$maximal_reads_number" '{arr[$1] += $3} END {for (key in arr) {if (arr[key]>0 && arr[key]<=var) print key}}' $out_demtx/fastq.non.zero.reads > $out_demtx/filtered.fastq.name # step 2 prune based on reads R1+R2
fi

echo -e "\n removing samples based on reads number finished"

################################################################################################################
############################################## quality control using trimmomatic  ##############################
echo -e "\nQuality control\n"

mkdir -p $parameter_o/trimmomatic

sed -e "s|^|$out_demtx/barcode3/|" -e "s|$|_R1.fastq|" $out_demtx/filtered.fastq.name|while read -r f1;\
do f2_file=$(basename $f1|cut -d_ -f1)_R2.fastq;dir_path=$(dirname $f1);f2="$dir_path/$f2_file";sample=$(basename $f1|cut -d_ -f1);\
$trimmomatic_path PE -threads $parameter_t $f1 $f2 $parameter_o/trimmomatic/${sample}_R1.paired.fastq $parameter_o/trimmomatic/${sample}_R1.unpaired.fastq $parameter_o/trimmomatic/${sample}_R2.paired.fastq $parameter_o/trimmomatic/${sample}_R2.unpaired.fastq \
ILLUMINACLIP:$parameter_primer:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:$parameter_trimmomatic_quality MINLEN:30 >> $parameter_o/trimmomatic/trimmomatic.log 2>&1;\
done
rm -f $parameter_o/filtered.fastq.name
echo -e "\nQuality control finished\n"
################################################################################################################
############################################## mapping reads to reference  #####################################
echo -e "\nmapping reads to reference\n"
mkdir -p $parameter_o/bwa_mem/index
ls $parameter_r/*|sort|while read -r line;do sample=$(basename $line|cut -d_ -f1); $bwa_mem2_path index -p $parameter_o/bwa_mem/index/${sample} $line;done

# method 1
mkdir -p $parameter_o/bwa_mem/sam
process_file2() {
  input=$1
  sample_fastq=$(basename "$input" | cut -d_ -f1 | cut -d- -f1)
  ref=$(find "$parallel_parameter_o/bwa_mem/index/" -name "*$sample_fastq*" | head -n 1)
  sample_barcode=$(basename "$input" | cut -d_ -f1)
  fastq1=$input
  fastq2=$(echo $fastq1 | sed "s/R1.paired.fastq/R2.paired.fastq/");
  if [[ -f "$ref" ]]; then
      $parallel_bwa_mem2_path mem -o $parallel_parameter_o/bwa_mem/sam/${sample_barcode}.sam -t 1 $parallel_parameter_o/bwa_mem/index/$sample_fastq $fastq1 $fastq2 >> $parallel_parameter_o/bwa_mem/sam/${sample_fastq}.log 2>&1
  else
      echo "Reference genome for $sample_fastq not found." >&2
  fi
}
export -f process_file2
export parallel_parameter_o="$parameter_o" # variable within parallel has to be exported
export parallel_bwa_mem2_path="$bwa_mem2_path"
time find "$parameter_o/trimmomatic/" -name '*R1.paired.fastq' | parallel -j $parameter_t -k process_file2 # parallel with 1 cpu is faster than non-parallel with multiple cpu 
echo -e "\nmapping reads to reference finished\n"

################################################################################################################
############################################## convert sam to bam, sort bam  ###################################
echo -e "\nconvert sam to bam\n"
mkdir -p $parameter_o/bwa_mem/bam
ls $parameter_o/bwa_mem/sam/*sam | parallel -j $parameter_t "sample=\$(basename {} | awk -F'.sam' '{print \$1}'); $samtools_path view --threads 1 -Sb {} > $parameter_o/bwa_mem/bam/\${sample}.bam; $samtools_path sort --threads 1 -o $parameter_o/bwa_mem/bam/\${sample}.sorted.bam $parameter_o/bwa_mem/bam/\${sample}.bam"
echo -e "\nconvert sam to bam finished\n"
################################################################################################################
###################################### calculate coverage for each position  ###################################
echo -e "calculate coverage for each position\n"
mkdir -p $parameter_o/coverage
ls $parameter_o/bwa_mem/bam/*sorted.bam|parallel -j $parameter_t " sample=\$(basename {}|awk -F".sorted.bam" '{print \$1}'); $samtools_path depth -a {} > $parameter_o/coverage/\${sample}.coverage.txt"
echo -e "calculate coverage for each position finished\n"

################################################################################################################
###################################### depth calculation pooling all barcodes  ###################################
echo -e "\n calculating pooled barcode position depth\n"

# summing up reads from the same contigs and same position , combined depth for each position 
awk '{data[$1" "$2] += $3}END{for (key in data) print key,data[key]}' $parameter_o/coverage/*.coverage.txt > $parameter_o/coverage/combined.depth.txt

# combined coverage
awk '{if ($3>0) data[$1] = data[$1] + 1; row[$1] = row[$1] + 1} END {for (key in data) {cov = (data[key]/row[key]);print key, cov}}' $parameter_o/coverage/combined.depth.txt > $parameter_o/coverage/combined.coverage.txt

################################################################################################################
###################################### average depth and coverage for each barcode  ############################
echo -e "\n calculate average coverage for each barcode"
count_coverage() {
    input=$1
    cov=$(awk '{if ($3>0) sum = sum + 1}END{print sum/NR}' $input)
    echo $input $cov >> $parallel_parameter_o/coverage/average.coverage.txt
}
export -f count_coverage
export parallel_parameter_o=$parameter_o

ls -lhS $parameter_o/coverage/*.coverage.txt | grep -v 'combined'|grep -v 'average'|awk '$5 != 0 {print $9}' > $parameter_o/coverage/fastq.non.zero

parallel -j $parameter_t count_coverage :::: $parameter_o/coverage/fastq.non.zero

sed 's/.*\///' $parallel_parameter_o/coverage/average.coverage.txt | awk '{split($1, arr, "-"); sum[arr[1]] += $2; count[arr[1]] += 1} END {for (key in sum) print key, sum[key] / count[key]}' >> $parallel_parameter_o/coverage/mean.average.coverage.txt

echo -e "\n calculate average coverage for each barcode finished "
