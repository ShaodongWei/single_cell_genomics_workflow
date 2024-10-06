my_cutadapt() {
    # Check if cutadapt is installed
    if ! command -v cutadapt &> /dev/null; then
        echo "Error: cutadapt is not installed or not available in the PATH."
        return 1
    fi

    # Initialize default values
    fastq1=""
    fastq2=""
    barcode=""
    output_dir=""
    threads=""
    error_rate=""
    overlap_minlength=""

    # Parse labeled parameters
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --fastq1) fastq1="$2"; shift ;;
            --fastq2) fastq2="$2"; shift ;;
            --barcode) barcode="$2"; shift ;;
            --output_dir) output_dir="$2"; shift ;;
            --threads) threads="$2"; shift ;;
            --error_rate) error_rate="$2"; shift ;;
            --overlap_minlength) overlap_minlength="$2"; shift ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if required arguments are provided
    if [[ -z "$fastq1" || -z "$barcode" || -z "$output_dir" || -z "$threads" ]]; then
        echo "Usage: my_cutadapt --fastq1 <fastq1> --fastq2 <fastq2> --barcode <barcode> --output_dir <output_dir> --threads <threads> --error_rate <error_rate> --overlap_minlength <overlap_minlength>"
        return 1
    fi

    # check if directory exists, if not, create it 
    if [[ ! -d $output_dir ]]; then 
        mkdir -p $output_dir
    fi 

    sample=$(basename $fastq1 | cut -d_ -f1)
    barcode_name=$(basename $barcode | cut -d. -f1)
    
    cutadapt --action trim -e $error_rate --overlap $overlap_minlength --cores $threads \
        -g file:$barcode \
        -o $output_dir/"${sample}-{name}_R1.fastq" \
        -p $output_dir/"${sample}-{name}_R2.fastq" \
        $fastq2 $fastq1 >> $output_dir/${barcode_name}.log 2>&1
}


parallel_cutadapt(){
    # Define a list of software to check
    required_software=("cutadapt")

    # Check if each software is installed
    for software in "${required_software[@]}"; do
        if ! command -v $software &> /dev/null; then
            echo "Error: $software is not installed or not available in the PATH."
            return 1
        fi
    done

    # Initialize default values
    fastq1=""
    barcode=""
    output_dir=""
    threads=""
    error_rate=""
    overlap_minlength=""

    # Parse labeled parameters
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --fastq1) fastq1="$2"; shift ;;
            --barcode) barcode="$2"; shift ;;
            --output_dir) output_dir="$2"; shift ;;
            --threads) threads="$2"; shift ;;
            --error_rate) error_rate="$2"; shift ;;
            --overlap_minlength) overlap_minlength="$2"; shift ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if required arguments are provided
    if [[ -z "$fastq1" || -z "$barcode" || -z "$output_dir" || -z "$threads" ]]; then
        echo "Usage: parallel_cutadapt --fastq1 <fastq1> --barcode <barcode> --output_dir <output_dir> --threads <threads> --error_rate <error_rate> --overlap_minlength <overlap_minlength>"
        return 1
    fi

    # check if directory exists, if not, create it 
    if [[ ! -d $output_dir ]]; then 
        mkdir -p $output_dir
    fi 

    fastq2=$(echo $fastq1|sed 's|R1.fastq|R2.fastq|')
    sample=$(basename $fastq1|cut -d_ -f1)
    barcode_name=$(basename $barcode|cut -d. -f1)
    cutadapt --action trim -e $error_rate --overlap 8 --cores $threads \
        -g file:$barcode \
        -o $output_dir/${sample}+{name}_R2.fastq \
        -p $output_dir/${sample}+{name}_R1.fastq \
        $fastq2 $fastq1 >> $output_dir/${barcode_name}.log 2>&1
    }


count_reads(){
    # Define a list of software to check
    required_software=("bioawk")

    # Check if each software is installed
    for software in "${required_software[@]}"; do
        if ! command -v $software &> /dev/null; then
            echo "Error: $software is not installed or not available in the PATH."
            return 1
        fi
    done

    # Initialize default values
    fastx=""

    # Parse labeled parameters
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --fastx) fastx="$2"; shift ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if required arguments are provided
    if [[ -z "$fastx" ]]; then
        echo "Usage: count_reads --fastx <fastx> "
        return 1
    fi

    reads_num=$(bioawk -cfastx '{if (length($seq)>0) sum=sum+1}END{print sum}' $fastx) 
    if [[ "$reads_num" -gt 0 ]]; then
        output=$(basename "$fastx")
        echo $output $reads_num
    fi

}

prune_sample(){
    # input is a dataframe containing 3 columns, 1st column is the sample name, 2nd column is the file name, 3rd column is the reads number
    # the reads here is sum of R1 and R2 file in corresponding samples. 

    # Define a list of software to check
    required_software=("awk")

    # Check if each software is installed
    for software in "${required_software[@]}"; do
        if ! command -v $software &> /dev/null; then
            echo "Error: $software is not installed or not available in the PATH."
            return 1
        fi
    done

    # Initialize default values
    file=""
    min_reads=""
    max_reads=""
    

    # Parse labeled parameters
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --file) file="$2"; shift ;;
            --min_reads) min_reads="$2"; shift ;;
            --max_reads) max_reads="$2"; shift ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if required arguments are provided
    if [[ -z "$file" || (-z "$min_reads" && -z "$max_reads") ]]; then
        echo "Usage: prune_sample --file <file> --min_reads <min_reads> --max_reads <max_reads>"
        return 1
    fi

    if [[ ! -z "$min_reads" && ! -z "$max_reads" ]];then 
        awk -v var1="$min_reads" -v var2="$max_reads" '{arr[$1] += $3} END {for (key in arr) {if (arr[key]>=var1 && arr[key]<=var2) print key}}' $file # between min and max 
    elif [[ ! -z "$min_reads" && -z "$max_reads" ]]; then 
        awk -v var="$min_reads" '{arr[$1] += $3} END {for (key in arr) {if (arr[key]>=var) print key}}' $file # larger than min
    elif [[ -z "$min_reads" && ! -z "$max_reads" ]]; then 
        awk -v var="$max_reads" '{arr[$1] += $3} END {for (key in arr) {if (arr[key]>0 && arr[key]<=var) print key}}' $file # less than max
    fi
}

bwa_mapping(){
    # Define a list of software to check
    required_software=("bwa-mem2")

    # Check if each software is installed
    for software in "${required_software[@]}"; do
        if ! command -v $software &> /dev/null; then
            echo "Error: $software is not installed or not available in the PATH."
            return 1
        fi
    done

    # Initialize default values
    fastq1=""
    reference_index=""
    output_dir=""
    threads=""

    # Parse labeled parameters
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --fastq1) fastq1="$2"; shift ;;
            --reference_index) reference_index="$2"; shift ;;
            --output_dir) output_dir="$2"; shift ;;
            --threads) threads="$2"; shift ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if required arguments are provided
    if [[ -z "$fastq1" || -z "$reference_index" || -z "$output_dir" || -z "$threads" ]]; then
        echo "Usage: bwa_mapping --fastq1 <fastq1> --reference_index <reference_index> --output_dir <output_dir> --threads <threads>"
        return 1
    fi

    # if specified directory does not exist, then create it 
    if [[ ! -d $output_dir ]];then mkdir -p $output_dir;fi 

    sample_fastq=$(basename "$fastq1" | cut -d_ -f1 | cut -d- -f1)
    sample_barcode=$(basename "$fastq1" | cut -d_ -f1)
    fastq2=$(echo $fastq1 | sed "s/R1.paired.fastq/R2.paired.fastq/")
    bwa-mem2 mem -o $output_dir/${sample_barcode}.sam -t $threads $reference_index $fastq1 $fastq2 
    
    # if the mapping not successful for samples, e.g. due to meomory problem somehow (e.g. core dumped files)
    if [ $? -ne 0 ]; then
        echo "$(date): bwa-mem2 mem failed for ${sample_barcode}" >> $output_dir/error.log
        # Optionally: take further action like deleting corrupted outputs
        rm -f $output_dir/${sample_barcode}.sam
    fi
}

kma_mapping(){
    # Define a list of software to check
    required_software=("kma")

    # Check if each software is installed
    for software in "${required_software[@]}"; do
        if ! command -v $software &> /dev/null; then
            echo "Error: $software is not installed or not available in the PATH."
            return 1
        fi
    done

    # Initialize default values
    fastq1=""
    reference_index=""
    output_dir=""
    threads=""

    # Parse labeled parameters
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --fastq1) fastq1="$2"; shift ;;
            --reference_index) reference_index="$2"; shift ;;
            --output_dir) output_dir="$2"; shift ;;
            --threads) threads="$2"; shift ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if required arguments are provided
    if [[ -z "$fastq1" || -z "$reference_index" || -z "$threads" || -z "$output_dir" ]]; then
        echo "Usage: kma_mapping --fastq1 <fastq1> --reference_index <reference_index> --output_dir <output_dir> --threads <threads>"
        return 1
    fi

    sample_fastq=$(basename "$fastq1" | cut -d_ -f1 | cut -d- -f1)
    sample_barcode=$(basename "$fastq1" | cut -d_ -f1)
    fastq2=$(echo $fastq1 | sed "s/R1.paired.fastq/R2.paired.fastq/")
    kma -ipe $fastq1 $fastq2 -t_db $reference_index -o $output_dir -t $threads -sam -k 10 -nc -nf
    
    # if the mapping not successful for samples, e.g. due to meomory problem somehow (e.g. core dumped files)
    if [ $? -ne 0 ]; then
        echo "$(date): kma mapping failed for ${sample_barcode}" >> $output_dir/error.log
        # Optionally: take further action like deleting corrupted outputs
        #rm -f $output_dir/${sample_barcode}.sam
    fi
}
count_coverage(){
    # Define a list of software to check
    required_software=("bioawk")

    # Check if each software is installed
    for software in "${required_software[@]}"; do
        if ! command -v $software &> /dev/null; then
            echo "Error: $software is not installed or not available in the PATH."
            return 1
        fi
    done

    # Initialize default values
    input=""
    reference=""
    
    # Parse labeled parameters
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --input) input="$2"; shift ;;
            --reference) reference="$2"; shift ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if required arguments are provided
    if [[ -z "$input" || -z "$reference" ]]; then
        echo "Usage: count_coverage --input <input> --reference <reference>"
        return 1
    fi

    # calculation (can handle genome (single length) or plasmid (multiple length))
    cat $reference | bioawk -cfastx '{print $name"~"length($seq)}' | while read -r line;do
        sample=$(basename $input)
        sample=${sample%.depth*}
        ref_name=$(echo $line|cut -d'~' -f1)
        ref_length=$(echo $line|cut -d'~' -f2)
        # Get the number of lines in the input file for the corresponding reference (important for plasmid)
        row=$(cat $input | grep "$ref_name" | wc -l)
        echo $sample $ref_name $ref_length $row|awk '{print $1,$2,$4/$3}'
    done
    
}
