
my_cutadapt() {
    # Check if cutadapt is installed
    if ! command -v cutadapt &> /dev/null; then
        echo "Error: cutadapt is not installed or not available in the PATH."
        return 1
    fi

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
        echo "Usage: my_cutadapt --fastq1 <fastq1> --barcode <barcode> --output_dir <output_dir> --threads <threads> --error_rate <error_rate> --overlap_minlength <overlap_minlength>"
        return 1
    fi

    # check if directory exists, if not, create it 
    if [[ ! -d $output_dir ]]; then 
        mkdir -p $output_dir
    fi 

    fastq2=$(echo $fastq1 | sed 's|R1.fastq|R2.fastq|')
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

mapping(){
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
        echo "Usage: mapping --fastq1 <fastq1> --reference_index <reference_index> --output_dir <output_dir> --threads <threads>"
        return 1
    fi

    # if specified directory does not exist, then create it 
    if [[ ! -d $output_dir ]];then mkdir -p $output_dir;fi 

    sample_fastq=$(basename "$fastq1" | cut -d_ -f1 | cut -d- -f1)
    sample_barcode=$(basename "$fastq1" | cut -d_ -f1)
    fastq2=$(echo $fastq1 | sed "s/R1.paired.fastq/R2.paired.fastq/")
    bwa-mem2 mem -o $output_dir/${sample_barcode}.sam -t $threads $reference_index $fastq1 $fastq2 
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
    fastq=""
    reference=""
    
    # Parse labeled parameters
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --fastq) fastq="$2"; shift ;;
            --reference) reference="$2"; shift ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if required arguments are provided
    if [[ -z "$fastq" || -z "$reference" ]]; then
        echo "Usage: count_coverage --fastq <fastq> --reference <reference>"
        return 1
    fi

    input=$1
    
    # Get reference length
    ref_len=$(bioawk -cfastx '{print length($seq)}' "$reference")

    # Get the number of lines in the input file
    row=$(wc -l < "$fastq" | awk '{print $1}')

    # Calculate coverage
    cov=$(awk -v row="$row" -v len="$ref_len" 'BEGIN { print row / len }')

    # Append coverage information to the output file
    echo "$fastq $cov"
}
  
