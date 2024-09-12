#!/bin/bash


# problem: so far the pipeline can only output results for each sample separately in separately output folder.

# version 2.13: (1) changed naming to depth for each barcode, 
               
################################################################################################################
############################################## start counting time #############################################
start=$(date +%s)


################################################################################################################
############################################## load modules from computerokme ##################################
module load tools computerome_utils/2.0 #load mudles from Computerome 2.0
module load parallel/20220422
module load samtools/1.18
module load bioawk/1.0
module load bwa-mem2/2.2.1
module load jdk/22
module load trimmomatic/0.38
module load anaconda3/2023.09-0

################################################################################################################
############################################## help page #######################################################
Help()
{
   # Display Help
   echo -e "backup your files first !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
   echo "Add your parameters here."
   echo "options:"
   echo "-i|--input   Necessary. Path to the raw undemultiplexed fastq/fastq.gz files directory (only one pair of fastq accepted). Sample names have to be named as by 'sample_R1*' or sample_R2*."
   echo "-r|--reference  Necessary. Path to the reference file path. " 
   echo "-o|--output    Necessary. Path to the output directory"
   echo "--barcode1 Necessary. The barcode1.fasta file path. The barcode1 is the one closest to the reads end."
   echo "--barcode2 Necessary. The barcode2.fasta file path. The barcode2 is the one 2nd closest to the reads end."
   echo "--barcode3 Necessary. The barcode3.fasta file path. The barcode3 is the one 3rd closest to the reads end."
   echo "--primer  Necessary. The primer fasta file path."
   echo "--minimal_reads_number   Optional. The minimal number of reads R1+R2, a value or 'None'. "
   echo "--maximal_reads_number   Optional. The maximal number of reads R1+R2, a value or 'None'. "
   echo "-t|--threads   Necessary. Number of threads to use"
   echo "--trimmomatic_quality  Optional. The quality score for trimmomatic windows."
   echo "--memory   If use "--memory", you will specify the memory to be used (GB)"

}

################################################################################################################
############################################## set parameters (for developer) ##################################

# initialize parameters
parameter_memory=20 #in GB
parameter_t=40
parameter_trimmomatic_quality=15

##
VALID_ARGS=$(getopt -o hi:r:o:t: --long input:reference:,output:,threads:,barcode1:,barcode2:,barcode3:,primer:,memory:,\
trimmomatic_quality:,gzip:,minimal_reads_number:,maximal_reads_number: -- "$@")
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
if [[ ! -z $parameter_r && -f $parameter_r ]]; then echo -e "\nYour reference for mapping is $parameter_r"; else echo -e "\nYour reference for mapping is missing !!!"; echo "use -h for help"; exit;fi
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

echo "You will use $parameter_t threads"
echo "You will use $parameter_memory GB memory"
echo "You will use trimmomatic quality $parameter_trimmomatic_quality"

################################################################################################################
############################################## file renaming  ##################################################
first_file=$(ls $parameter_i|head -n 1)
check=$(echo $first_file|grep -e '_R1.fastq' -e '_R2.fastq')
if [[ -z $check ]];then
  echo -e "\n Renaming files\n"
  ls $parameter_i|while read -r line;do
  new_name1=$(echo $line|sed -e 's/\.R1.fastq/_R1.fastq/' -e 's/\.R1.fq/_R1.fastq/' -e 's/\.1.fastq/_R1.fastq/' -e 's/\.1.fq/_R1.fastq/');\
  new_name2=$(echo $new_name1|sed -e 's/\.R2.fastq/_R2.fastq/' -e 's/\.R2.fq/_R2.fastq/' -e 's/\.2.fastq/_R2.fastq/' -e 's/\.2.fq/_R2.fastq/');\
  if [[ $line != $new_name2 ]];then
  mv $parameter_i/$line $parameter_i/$new_name2;\
  fi;
  done
  echo -e "\n Renaming files finished\n"
fi


################################################################################################################
############################################## demultiplexing  ##################################################
# check if demultiplexed dir already exists
if [[ ! -d $parameter_o/demultiplexed ]]; then
  echo -e "\nDemultiplexing\n"
  mkdir $parameter_o/demultiplexed
  mkdir $parameter_o/demultiplexed/barcode1
  mkdir $parameter_o/demultiplexed/barcode2
  mkdir $parameter_o/demultiplexed/barcode3

  # I think the first demultiplexing step should use all cores, from 2nd step, since we have >40 barcode, so we then use 1 core
  parallel_cutadapt1(){
    fastq1=$1
    fastq2=$(echo $fastq1|sed 's|R1.fastq|R2.fastq|')
    sample=$(basename $fastq1|cut -d_ -f1)
    barcode_name=$(basename $parallel_parameter_b1|cut -d. -f1)
    cutadapt --action trim -e 0.125 --overlap 8 --cores $parallel_threads -g file:$parallel_parameter_b1 -o $parallel_out/demultiplexed/barcode1/${sample}-{name}_R2.fastq -p $parallel_out/demultiplexed/barcode1/${sample}-{name}_R1.fastq $fastq2 $fastq1 >> $parallel_out/demultiplexed/barcode1/${barcode_name}.log 2>&1
    }
  export parallel_parameter_b1=$parameter_b1
  export parallel_out=$parameter_o
  export parallel_threads=$parameter_t
  export -f parallel_cutadapt1
  parallel -j 1 parallel_cutadapt1 ::: $parameter_i/*R1*

  parallel_cutadapt2(){
    fastq1=$1
    fastq2=$(echo $fastq1|sed 's|R1.fastq|R2.fastq|')
    sample=$(basename $fastq1|cut -d_ -f1)
    barcode_name=$(basename $parallel_parameter_b2|cut -d. -f1)
    cutadapt --action trim -e 0.125 --overlap 8 --cores 1 -g file:$parallel_parameter_b2 -o $parallel_out/demultiplexed/barcode2/${sample}+{name}_R2.fastq -p $parallel_out/demultiplexed/barcode2/${sample}+{name}_R1.fastq $fastq2 $fastq1 >> $parallel_out/demultiplexed/barcode2/${barcode_name}.log 2>&1
    }
  export parallel_parameter_b2=$parameter_b2
  export parallel_out=$parameter_o
  export -f parallel_cutadapt2
  find $parameter_o/demultiplexed/barcode1 -type f -name "*-*R1.fastq" -size +0 ! -name '*unknown*' > $parameter_o/demultiplexed/demultiplexed.barcode1.file.size.not.zero.fastq # to skip zero reads files
  #ls -lhS $parameter_o/demultiplexed/barcode1/*|grep 'R1.fastq'|grep -v 'unknown'|awk '$5 != 0 {print $9}' > $parameter_o/demultiplexed/demultiplexed.barcode1.file.size.not.zero.fastq # to skip zero reads files
  parallel -j $parameter_t parallel_cutadapt2 :::: $parameter_o/demultiplexed/demultiplexed.barcode1.file.size.not.zero.fastq
  #rm $parameter_o/demultiplexed/demultiplexed.barcode1.file.size.not.zero.fastq

  parallel_cutadapt3(){
    fastq1=$1
    fastq2=$(echo $fastq1|sed 's|R1.fastq|R2.fastq|')
    sample=$(basename $fastq1|cut -d_ -f1)
    barcode_name=$(basename $parallel_parameter_b3|cut -d. -f1)
    cutadapt --action trim -e 0.125 --overlap 8 --cores 1 -g file:$parallel_parameter_b3 -o $parallel_out/demultiplexed/barcode3/${sample}+{name}_R2.fastq -p $parallel_out/demultiplexed/barcode3/${sample}+{name}_R1.fastq $fastq2 $fastq1 >> $parallel_out/demultiplexed/barcode3/${barcode_name}.log 2>&1
    }
  export parallel_parameter_b3=$parameter_b3
  export parallel_out=$parameter_o
  export -f parallel_cutadapt3
  find $parameter_o/demultiplexed/barcode2 -type f -name "*-*R1.fastq" -size +0 ! -name '*unknown*' > $parameter_o/demultiplexed/demultiplexed.barcode2.file.size.not.zero.fastq # to skip zero reads files
  #ls -lhS $parameter_o/demultiplexed/barcode2/*|grep 'R1.fastq'| grep -v 'unknown'|awk '$5 != 0 {print $9}' > $parameter_o/demultiplexed/demultiplexed.barcode2.file.size.not.zero.fastq # to skip zero reads files
  parallel -j $parameter_t parallel_cutadapt3 :::: $parameter_o/demultiplexed/demultiplexed.barcode2.file.size.not.zero.fastq
  #rm $parameter_o/demultiplexed/demultiplexed.barcode2.file.size.not.zero.fastq
  echo -e "\nDemultiplexing finished\n"
fi

################################################################################################################
############################################## remove samples based on reads number  ###########################
# prune based on file size
if [[ ! -f $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads ]]; then
    echo -e "\nremoving samples based on reads number"
    find $parameter_o/demultiplexed/barcode3 -type f -name "*-*fastq" -size +0 ! -name '*unknown*' > $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq # prune based on file size
    #ls -lhS "$parameter_o/demultiplexed/barcode3"|grep -v 'unknown'|awk '$5 != 0 {print $9}'|sed '1d' > $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq # prune based on file size
    #sed -i '/barcode[1-3].log/d' $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq

    # calculate reads for each fastq file
    process_file() {
    input=$1
    reads_num=$(bioawk -cfastx '{if (length($seq)>0) sum=sum+1}END{print sum}' $input) #we count reads when length > 0, sometimes it can have reads name but not sequences so grep -c '^@' will always give same R1 R2 reads number, but our method not necessary R1 always equal R2 number.
    if [[ "$reads_num" -gt 0 ]]; then
    output=$(basename "$input")
    echo $output $reads_num >> $parallel_out/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads
    fi
    }
    export -f process_file # in process_file we keep all reads > 0
    export parallel_out="$parameter_o"
    #sed -i "s|^|$parameter_o/demultiplexed/barcode3/|" $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq # add prefix path
    parallel -j $parameter_t process_file :::: $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq
    awk '{split($1,arr,"_");print arr[1],$0}' $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads > $parameter_o/demultiplexed/tmp && mv $parameter_o/demultiplexed/tmp $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads
fi

# prune based on reads R1+R2
if [[ ! -f $parameter_o/demultiplexed/filtered.fastq.name ]]; then
  if [[ -n $minimal_reads_number && -n $maximal_reads_number ]]; then
          awk -v var1="$minimal_reads_number" -v var2="$maximal_reads_number" '{arr[$1] += $3} END {for (key in arr) {if (arr[key]>=var1 && arr[key]<=var2) print key}}' $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads > $parameter_o/demultiplexed/filtered.fastq.name # between min and max
      elif [[ -n $minimal_reads_number && -z $maximal_reads_number ]]; then
          awk -v var="$minimal_reads_number" '{arr[$1] += $3} END {for (key in arr) {if (arr[key]>=var) print key}}' $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads > $parameter_o/demultiplexed/filtered.fastq.name # larger than min
      elif [[ -n $maximal_reads_number && -z $minimal_reads_number ]]; then
          awk -v var="$maximal_reads_number" '{arr[$1] += $3} END {for (key in arr) {if (arr[key]>0 && arr[key]<=var) print key}}' $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads > $parameter_o/demultiplexed/filtered.fastq.name # less than max
  fi
echo -e "\nremoving samples based on reads number finished"
fi

################################################################################################################
############################################## calculate percentage of reads demultiplexed  ####################

if [[ ! -f $parameter_o/demultiplexed/raw.and.demultiplexed.reads ]]; then
  echo -e "\ncalculating percentage of reads demultiplexed"

  # Function to process each file
  count_reads() {
      local file_path=$1
      local total_lines=$(wc -l < "$file_path")
      local total_reads=$((total_lines / 4))
      local filename=$(basename "$file_path")
      echo "$filename $total_reads"
  }

  export -f count_reads

  # Use GNU Parallel to process files in parallel
  find "$parameter_i" -name '*.fastq' | parallel count_reads {} >> "$parameter_o/demultiplexed/sample.raw.reads"

  # sum up demultiplexed reads for each sample
  awk '{split($1,a,"-"); split($2,b,"_"); print a[1]"_"b[2]"\t"$NF}' $parameter_o/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads \
  |awk '{arr[$1] += $2} END {for (key in arr) print key"\t"arr[key]}'| sort -k1 > $parameter_o/demultiplexed/demultiplexed.reads

  # calculate proportion of demultiplexed compared with raw reads.
  sort -k 1 $parameter_o/demultiplexed/sample.raw.reads > tmp && mv tmp $parameter_o/demultiplexed/sample.raw.reads
  join -1 1 -2 1 $parameter_o/demultiplexed/sample.raw.reads $parameter_o/demultiplexed/demultiplexed.reads | awk '{print $0,$3/$2}' > $parameter_o/demultiplexed/raw.and.demultiplexed.reads
  echo -e "\ncalculating percentage of reads demultiplexed finished"
fi


################################################################################################################
############################################## quality control using trimmomatic  ##############################
if [[ ! -d $parameter_o/trimmomatic ]]; then
  echo -e "\nQuality control using trimmomatic\n"
  mkdir $parameter_o/trimmomatic
  sed -e "s|^|$parameter_o/demultiplexed/barcode3/|" -e "s|$|_R1.fastq|" $parameter_o/demultiplexed/filtered.fastq.name|while read -r f1;\
  do f2_file=$(basename $f1|cut -d_ -f1)_R2.fastq;dir_path=$(dirname $f1);f2="$dir_path/$f2_file";sample=$(basename $f1|cut -d_ -f1);\
  trimmomatic PE -threads $parameter_t $f1 $f2 $parameter_o/trimmomatic/${sample}_R1.paired.fastq $parameter_o/trimmomatic/${sample}_R1.unpaired.fastq $parameter_o/trimmomatic/${sample}_R2.paired.fastq $parameter_o/trimmomatic/${sample}_R2.unpaired.fastq \
  ILLUMINACLIP:$parameter_primer:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:$parameter_trimmomatic_quality MINLEN:30 >> $parameter_o/trimmomatic/trimmomatic.log 2>&1;\
  done
  echo -e "\nQuality control finished\n"
fi

################################################################################################################
############################################## mapping reads to reference  #####################################

if [[ ! -d $parameter_o/bwa_mem ]];then
  echo -e "\nmapping reads to reference\n"
  mkdir -p $parameter_o/bwa_mem/index

  # index the reference file 
  sample=$(basename "$parameter_r" | cut -d_ -f1)
  bwa-mem2 index -p $parameter_o/bwa_mem/index/$sample $parameter_r > /dev/null


  # method 1
  mkdir -p $parameter_o/bwa_mem/sam
  process_file2() {
    input=$1
    sample_fastq=$(basename "$input" | cut -d_ -f1 | cut -d- -f1)
    ref=$parallel_parameter_r
    sample_barcode=$(basename "$input" | cut -d_ -f1)
    fastq1=$input
    fastq2=$(echo $fastq1 | sed "s/R1.paired.fastq/R2.paired.fastq/");
    if [[ -f "$ref" ]]; then
        bwa-mem2 mem -o $parallel_parameter_o/bwa_mem/sam/${sample_barcode}.sam -t 1 $parallel_parameter_o/bwa_mem/index/$sample_fastq $fastq1 $fastq2 >> $parallel_parameter_o/bwa_mem/sam/${sample_fastq}.log 2>&1
    else
        echo "Reference genome for $sample_fastq not found." >&2
    fi
  }
  export -f process_file2
  export parallel_parameter_o="$parameter_o" # variable within parallel has to be exported
  export parallel_parameter_r="$parameter_r"
  # paired file is non-zero file size
  find "$parameter_o/trimmomatic/" -name '*R1.paired.fastq' | parallel -j $parameter_t -k process_file2 # parallel with 1 cpu is faster than non-parallel with multiple cpu
  echo -e "\nmapping reads to reference finished\n"
fi


################################################################################################################
############################################## convert sam to bam, sort bam  ###################################

if [[ ! -d $parameter_o/bwa_mem/bam ]]; then
  echo -e "\nconvert sam to bam\n"
  mkdir -p $parameter_o/bwa_mem/bam
  find $parameter_o/bwa_mem/sam -name *sam -type f | parallel -j $parameter_t "sample=\$(basename {} | awk -F'.sam' '{print \$1}'); samtools view --threads 1 -Sb {} > $parameter_o/bwa_mem/bam/\${sample}.bam; samtools sort --threads 1 -o $parameter_o/bwa_mem/bam/\${sample}.sorted.bam $parameter_o/bwa_mem/bam/\${sample}.bam"
  echo -e "\nconvert sam to bam finished\n"
fi
################################################################################################################
###################################### calculate depth for each barocde  #######################################
if [[ ! -d $parameter_o/coverage ]]; then
  echo -e "calculate coverage for each position\n"
  mkdir -p $parameter_o/coverage
  find $parameter_o/bwa_mem/bam -name *sorted.bam -type f | parallel -j $parameter_t " sample=\$(basename {}|awk -F".sorted.bam" '{print \$1}'); samtools depth {} > $parameter_o/coverage/\${sample}.depth.txt"
  echo -e "calculate depth for each position finished\n"
fi

################################################################################################################
############### calculate depth && coverage for each contig  ###################################################
if [[ ! -f $parameter_o/coverage/contigs.depth.each.position.txt ]]; then

  echo -e "\ncalculate depth for each position for each barcode\n"
  # summing up reads from the same contigs and same position , combined depth for each position from each barcode file
  find $parameter_o/coverage/ -type f -name "*.depth.txt" -print0 | xargs -0 cat | awk '{data[$1" "$2] += $3} END {for (key in data) print key, data[key]}' > $parameter_o/coverage/contigs.depth.each.position.txt

  # Calculate the length of reference and store it in a variable 
  ref_length=$( bioawk -cfastx '{print $name"\t"length($seq)}' "$parameter_r")
  # calculate the percentage of positions having reads mapped (genome coverage)
  awk -v ref_length="$ref_length" 'BEGIN{split(ref_length, parts, "\t"); total[parts[1]] = parts[2]} {if ($3 > 0) {data[$1] += 1}} END {for (key in data) {print key, data[key] / total[key]}}' $parameter_o/coverage/contigs.depth.each.position.txt > $parameter_o/coverage/contigs.coverage.txt
  
fi

################################################################################################################
###################################### coverage for each barcode && coverage for each sample  ##################
if [[ ! -f $parameter_o/coverage/fastq.file.size.not.zero ]]; then
  echo -e "\ncalculate coverage for each barcode"

  # calculate percentage of positions having reads mapped for each barcode (coverage)
  count_coverage() {
  input=$1
  # # Extract unique references (assuming this step is necessary)
  # ref=$(awk '{print $1}' "$input" | uniq)
  
  # Get reference length
  ref_len=$(bioawk -cfastx '{print length($seq)}' "$function_parameter_r")

  # Get the number of lines in the input file
  row=$(wc -l < "$input" | awk '{print $1}')

  # Calculate coverage
  cov=$(awk -v row="$row" -v len="$ref_len" 'BEGIN { print row / len }')

  # Append coverage information to the output file
  echo "$input $cov" >> "$function_parameter_o/coverage/barcode.coverage.txt"
  }
  export -f count_coverage
  export function_parameter_r="$parameter_r"
  export function_parameter_o="$parameter_o"

  find $parameter_o/coverage/ -type f -name "*-*.depth.txt" -size +0 > $parameter_o/coverage/fastq.file.size.not.zero
  parallel -j $parameter_t count_coverage :::: $parameter_o/coverage/fastq.file.size.not.zero

  # calculate coverage for each sample (mean of coverage from a sample)
  sed 's/.*\///' $parallel_parameter_o/coverage/barcode.coverage.txt | awk '{split($1, arr, "-"); sum[arr[1]] += $2; count[arr[1]] += 1} END {for (key in sum) print key, sum[key] / count[key]}' >> $parallel_parameter_o/coverage/sample.coverage.txt
  echo -e "\ncalculate average coverage for each barcode finished "
fi


################################################################################################################
############################################## finish counting time ################################################
end=$(date +%s)
runtime=$((end-start))
runtime_min=$(echo "scale=2; $runtime / 60" | bc)
echo "Script ran in $runtime_min minutes."





