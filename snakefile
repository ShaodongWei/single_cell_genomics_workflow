configfile: "config.yaml" 

rule all:
    input:
        config["output_directory"] + "/demultiplexed/.demultiplexing_barcode1.done",
        config["output_directory"] + "/demultiplexed/.demultiplexing_barcode2.done",
        config["output_directory"] + "/demultiplexed/.demultiplexing_barcode3.done",
        config['output_directory'] + '/demultiplexed/.prune_sample_by_reads.done',
        config['output_directory'] + '/demultiplexed/.reads_percentage_demultiplexed.done',
        config['output_directory'] + '/trimmomatic/.quality_control.done',
        config['output_directory'] + '/bwa_mem/.mapping2reference.done',
        config['output_directory'] + '/depth_coverage/.depth_coverage.done'
        
rule demultiplexing:
    input:
        fastq1 = config['fastq1'],
        fastq2 = config['fastq2']
    output:
        [config["output_directory"] + "/demultiplexed/.demultiplexing_barcode1.done",
        config["output_directory"] + "/demultiplexed/.demultiplexing_barcode2.done",
        config["output_directory"] + "/demultiplexed/.demultiplexing_barcode3.done"]
    params:
        output_dir = config['output_directory'],
        barcode1 = config["barcode3"],  
        barcode2 = config["barcode2"],
        barcode3 = config["barcode1"],
        threads = config['threads'],
        error_rate = config['cutadapt_error_rate'],
        overlap_minlength = config['cutadapt_overlap_minlength'],
        functions = config["functions"]
    conda:
        "env/demultiplexing.yaml"
    shell:
        '''
        echo {params.output_dir}
        mkdir -p {params.output_dir}/demultiplexed
        mkdir -p {params.output_dir}/demultiplexed/barcode1
        mkdir -p {params.output_dir}/demultiplexed/barcode2
        mkdir -p {params.output_dir}/demultiplexed/barcode3

        # load functions 
        source {params.functions} # load functions 

        # step1, the first demultiplexing step should use all cores, from 2nd step, since we have >40 barcode, so we then use 1 core in each parallel 
        
        my_cutadapt --fastq1 {input.fastq1} --fastq2 {input.fastq2} --barcode {params.barcode1} --output_dir {params.output_dir}/demultiplexed/barcode1 --threads {params.threads} --error_rate {params.error_rate} --overlap_minlength {params.overlap_minlength}
        
        touch {output[0]} # Flag to indicate completion, within curly bracket is python style 

        # step 2, demultiplex barcode 2 in parallel, that in each parallel using only 1 thread, 
        find {params.output_dir}/demultiplexed/barcode1 -type f -name "*-*R1.fastq" -size +0 ! -name '*unknown*' > {params.output_dir}/demultiplexed/demultiplexed.barcode1.file.size.not.zero.fastq # to skip zero reads files
        export -f parallel_cutadapt
        export TMPDIR=/tmp
        parallel -j {params.threads} "parallel_cutadapt --fastq1 {{}} --barcode {params.barcode2} --output_dir {params.output_dir}/demultiplexed/barcode2 --threads 1 --error_rate {params.error_rate} --overlap_minlength {params.overlap_minlength}" :::: {params.output_dir}/demultiplexed/demultiplexed.barcode1.file.size.not.zero.fastq
        touch {output[1]}

        # step 3, demultiplex barcode 3 in parallel, that in each parallel using only thread, 
        find {params.output_dir}/demultiplexed/barcode2 -type f -name "*-*R1.fastq" -size +0 ! -name '*unknown*' > {params.output_dir}/demultiplexed/demultiplexed.barcode2.file.size.not.zero.fastq # to skip zero reads files
        export -f parallel_cutadapt
        parallel -j {params.threads} "parallel_cutadapt --fastq1 {{}} --barcode {params.barcode3} --output_dir {params.output_dir}/demultiplexed/barcode3 --threads 1 --error_rate {params.error_rate} --overlap_minlength {params.overlap_minlength}" :::: {params.output_dir}/demultiplexed/demultiplexed.barcode2.file.size.not.zero.fastq
        
        touch {output[2]}

        '''

rule prune_sample_by_reads:
    input:
        config["output_directory"] + "/demultiplexed/.demultiplexing_barcode3.done"
    output:
        config['output_directory'] + '/demultiplexed/.prune_sample_by_reads.done'
    params:
        output_dir = config['output_directory'],
        threads = config['threads'],
        min_reads = config['min_reads'],
        max_reads = config['max_reads'],
        functions = config["functions"]
    conda:
        "env/prune_sample_by_reads.yaml"
    shell:
        '''
        source {params.functions}

        # select non zero fastq files
        find {params.output_dir}/demultiplexed/barcode3 -type f -name "*-*fastq" -size +0 ! -name '*unknown*' > {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq
        
        # calculate reads for each fastq file 
        export -f count_reads
        if [[ -f {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads ]];then 
            rm {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads;
        fi
        export TMPDIR=/tmp
        parallel -j {params.threads} "count_reads --fastx {{}} >> {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads" \
            :::: {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq
        
        # make a new column 'sample'
        awk '{{split($1,arr,"_");print arr[1],$0}}' {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads > {params.output_dir}/demultiplexed/tmp \
            && mv {params.output_dir}/demultiplexed/tmp {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads
        
        # prune sample based on R1+R2 reads
        prune_sample --file {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads --min_reads {params.min_reads} --max_reads {params.max_reads} > {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads.selected
        
        touch {output}
        '''

rule reads_percentage_demultiplexed:
    input:
        config['output_directory'] + '/demultiplexed/.prune_sample_by_reads.done'
    output:
        config['output_directory'] + '/demultiplexed/.reads_percentage_demultiplexed.done'
    params:
        fastq1 = config['fastq1'],
        fastq2 = config['fastq2'],
        output_dir = config['output_directory'],
        threads = config['threads'],
        functions = config["functions"]
    conda:
        "env/reads_percentage_demultiplexed.yaml"
    shell:
        '''
        source {params.functions}

        # count raw reads
        export -f count_reads
        if [[ -f {params.output_dir}/demultiplexed/sample.raw.reads ]]; then rm {params.output_dir}/demultiplexed/sample.raw.reads; fi
        export TMPDIR=/tmp # for parallel, temporary files, to avoid the system /etc /scratch where I have read-only permission. 
        echo -e "{params.fastq1}\n{params.fastq2}" | parallel -j {params.threads} "count_reads --fastx {{}} >> {params.output_dir}/demultiplexed/sample.raw.reads"

        # count demultiplexed reads for each file
        awk '{{split($1,a,"-"); split($2,b,"_"); print a[1]"_"b[2]"\t"$NF}}' {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads \
        |awk '{{arr[$1] += $2}} END {{for (key in arr) print key"\t"arr[key]}}'| sort -k1 > {params.output_dir}/demultiplexed/demultiplexed.reads

        # calculate proportion of demultiplexed compared with raw reads
        sort -k 1 {params.output_dir}/demultiplexed/sample.raw.reads > tmp && mv tmp {params.output_dir}/demultiplexed/sample.raw.reads
        join -1 1 -2 1 {params.output_dir}/demultiplexed/sample.raw.reads {params.output_dir}/demultiplexed/demultiplexed.reads | awk '{{print $0,$3/$2}}' \
            > {params.output_dir}/demultiplexed/raw.percentage.demultiplexed.reads

        touch {output}
        '''

rule quality_control:
    input:
        config['output_directory'] + '/demultiplexed/.reads_percentage_demultiplexed.done'
    output:
        config['output_directory'] + '/trimmomatic/.quality_control.done'
    params:
        output_dir = config['output_directory'],
        threads = config['threads'],
        primer = config['primer'],
        trimmomatic_quality = config['trimmomatic_quality']
    conda:
        "env/quality_control.yaml" 
    shell:
        '''
        mkdir -p {params.output_dir}/trimmomatic

        if [[ -f {params.output_dir}/trimmomatic/trimmomatic.log ]]; then rm {params.output_dir}/trimmomatic/trimmomatic.log; fi
        sed -e "s|^|{params.output_dir}/demultiplexed/barcode3/|" -e "s|$|_R1.fastq|" {params.output_dir}/demultiplexed/demultiplexed.barcode3.file.size.not.zero.fastq.reads.selected |while read -r f1;\
        do f2_file=$(basename $f1|cut -d_ -f1)_R2.fastq;dir_path=$(dirname $f1);f2="$dir_path/$f2_file";sample=$(basename $f1|cut -d_ -f1);\
        trimmomatic PE -threads {params.threads} $f1 $f2 {params.output_dir}/trimmomatic/${{sample}}_R1.paired.fastq {params.output_dir}/trimmomatic/${{sample}}_R1.unpaired.fastq \
            {params.output_dir}/trimmomatic/${{sample}}_R2.paired.fastq {params.output_dir}/trimmomatic/${{sample}}_R2.unpaired.fastq \
        ILLUMINACLIP:{params.primer}:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:{params.trimmomatic_quality} MINLEN:30 >> {params.output_dir}/trimmomatic/trimmomatic.log 2>&1;\
        done

        touch {output}

        '''

rule mapping:
    input:
        config['output_directory'] + '/trimmomatic/.quality_control.done'
    output:
        config['output_directory'] + '/bwa_mem/.mapping2reference.done'
    params:
        output_dir = config['output_directory'],
        threads = config['threads'],
        reference = config['reference'],
        functions = config['functions'],
    conda:
        "env/mapping.yaml" 

    shell:
        r'''              
        source {params.functions}
        
        mkdir -p {params.output_dir}/bwa_mem/index
        mkdir -p {params.output_dir}/bwa_mem/sam

        # index the reference file 
        sample=$(basename {params.reference} | cut -d_ -f1)
        bwa-mem2 index -p {params.output_dir}/bwa_mem/index/$sample {params.reference} > /dev/null
    
        # mapping
        if [[ -f {params.output_dir}/bwa_mem/sam/bwa-mem2.log ]];then rm {params.output_dir}/bwa_mem/sam/bwa-mem2.log;fi
        if [[ -f {params.output_dir}/bwa_mem/sam/error.log ]];then rm {params.output_dir}/bwa_mem/sam/error.log;fi
        export -f mapping
        export TMPDIR=/tmp # for parallel, temporary files, to avoid the system /etc /scratch where I have read-only permission. 

        find "{params.output_dir}/trimmomatic/" -name '*R1.paired.fastq' | parallel -j {params.threads} -k \
            "mapping --fastq1 {{}} --reference_index {params.output_dir}/bwa_mem/index/$sample --output_dir {params.output_dir}/bwa_mem/sam --threads 1 \
                >> {params.output_dir}/bwa_mem/sam/bwa-mem2.log 2>&1"

        # converting sam to bam 
        mkdir -p {params.output_dir}/bwa_mem/bam
        export TMPDIR=/tmp
        # we need to use r to escape the backslash 
        find {params.output_dir}/bwa_mem/sam -name *sam -type f | parallel -j {params.threads} "sample=\$(basename {{}} | awk -F'.sam' '{{print \$1}}'); samtools view --threads 1 -Sb {{}} > {params.output_dir}/bwa_mem/bam/\${{sample}}.bam; samtools sort --threads 1 -o {params.output_dir}/bwa_mem/bam/\${{sample}}.sorted.bam {params.output_dir}/bwa_mem/bam/\${{sample}}.bam"
        
        touch {output}
        '''

rule depth_coverage:
    input:
        config['output_directory'] + '/bwa_mem/.mapping2reference.done'
    output:
        config['output_directory'] + '/depth_coverage/.depth_coverage.done'
    params:
        output_dir = config['output_directory'],
        threads = config['threads'],
        reference = config['reference'],
        functions = config['functions'],
    conda:
        "env/depth_coverage.yaml" 
    shell:
        r'''
        mkdir -p {params.output_dir}/depth_coverage

        # sequencing depth for non-zero position for each barcode 
        export TMPDIR=/tmp
        find {params.output_dir}/bwa_mem/bam -name *sorted.bam -type f | parallel -j {params.threads} " sample=\$(basename {{}}|awk -F".sorted.bam" '{{print \$1}}'); \
            samtools depth {{}} > {params.output_dir}/depth_coverage/\${{sample}}.depth.txt"
        
        # depth for each contig, summing up reads from the same contigs and same position , combine depth for each position from each barcode file
        find {params.output_dir}/depth_coverage/ -type f -name "*.depth.txt" -print0 | xargs -0 cat | awk '{{data[$1" "$2] += $3}} END {{for (key in data) print key, data[key]}}' \
            > {params.output_dir}/depth_coverage/contigs.depth.each.position.txt

        # coverage for each contig, Calculate the length of reference and store it in a variable 
        ref_length=$( bioawk -cfastx '{{print $name"\t"length($seq)}}' {params.reference}) #when for plasmid, it can be multiple lines 
        echo "$ref_length"
        
        # calculate the percentage of positions having reads mapped (genome coverage)
        if [[ -f {params.output_dir}/depth_coverage/contigs.coverage.txt ]];then rm {params.output_dir}/depth_coverage/contigs.coverage.txt;fi
        
        echo "$ref_length" | while read -r line;do
            awk -v ref_length_line="$line" 'BEGIN{{split(ref_length_line, parts, "\t"); total[parts[1]] = parts[2]}} {{if ($3 > 0) {{data[$1] += 1}}}} END \
                {{for (key in total) {{print key, data[key] / total[key]}}}}' {params.output_dir}/depth_coverage/contigs.depth.each.position.txt >> {params.output_dir}/depth_coverage/contigs.coverage.txt
        done
        
        touch {output}
        '''
