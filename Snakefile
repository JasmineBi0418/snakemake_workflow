from snakemake.utils import min_version
min_version("5.3.0")

wildcard_constraints:
    genome_name="[^_]+"

configfile: "config.yml"
sample_id = config["sample_id"]  # Changed variable name
genome_name = config["genome_name"]
# Description: Rule to identify differentially expressed genes and create visualizations.

rule identify_DE_genes_w_visualization:  
    input:
        expand("count_data/{genome_name}/{sample}_counts.txt", genome_name=genome_name, sample=sample_id),
        expand("count_data/{genome_name}/{sample}_counts.txt.summary", genome_name=genome_name, sample=sample_id),
        expand("salmon_quant/{sample}", sample=sample_id),

rule fastqc:
    input:
        r1="data/raw_internal/{sample}_1.fq.gz",
        r2="data/raw_internal/{sample}_2.fq.gz",
    output:
        html1="results/rawQC/{sample}_1_fastqc.html",
        zip1="results/rawQC/{sample}_1_fastqc.zip",
        html2="results/rawQC/{sample}_2_fastqc.html",
        zip2="results/rawQC/{sample}_2_fastqc.zip",
    threads: 6
    params: 
        out_dir="results/rawQC/"
    message: """--- Quality check of raw data with FastQC before trimming.---"""
    log:
        "logs/qc_{sample}.log"
    shell:
        '''
        mkdir -p results/rawQC
        fastqc -t {threads}  -o {params.out_dir} {input.r1}
        fastqc -t {threads}  -o {params.out_dir} {input.r2}
        '''

rule multiqc:
    """
    Aggregate all FastQC reports into a MultiQC report.
    """
    input: 
       html1="results/rawQC/{sample}_1_fastqc.html",
        zip1="results/rawQC/{sample}_1_fastqc.zip",
        html2="results/rawQC/{sample}_2_fastqc.html",
        zip2="results/rawQC/{sample}_2_fastqc.zip",
    output:
        "results/multiqc_results/multiqc_report.html"
    message: """--- Aggregate all FastQC reports into a MultiQC report.---"""
    log:
        "logs/multiqc.log"
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p results/multiqc_results
        # Run multiQC and keep the html report
        multiqc {input} --outdir results/multiqc_results --no-data-dir --no-ansi --force
        """
rule trimming:
    input:
        html1="results/rawQC/{sample}_1_fastqc.html",
        zip1="results/rawQC/{sample}_1_fastqc.zip",
        html2="results/rawQC/{sample}_2_fastqc.html",
        zip2="results/rawQC/{sample}_2_fastqc.zip",
        r1="data/raw_internal/{sample}_1.fq.gz",
        r2="data/raw_internal/{sample}_2.fq.gz",
    output:
        trimmed_r1="data/trimmed_data/trimmed_{sample}_1.fq.gz",
        trimmed_r2="data/trimmed_data/trimmed_{sample}_2.fq.gz",
        json="data/trimmed_data/trimmed_{sample}.json",
        html="data/trimmed_data/trimmed_{sample}.html",
    message: """--- Automated trimming raw data using fastp tool.---"""
    log:
        "logs/trimming_{sample}.log"
    shell:
        '''
        mkdir -p data/trimmed_data
        fastp -i {input.r1} -I {input.r2} -o {output.trimmed_r1} -O {output.trimmed_r2} -j {output.json} -h {output.html}
        '''

rule get_genome_fasta:
    message:
        """
        ---Retrieve the sequence in fasta format for a genome.---
        """
    output:
        genome_fasta= "data/raw_external/{genome_name}.genome.fa.gz"
    params:
        fasta_path = lambda wildcards: config["genomes"][wildcards.genome_name]["fasta"]
    log:
        "logs/{genome_name}_genome_download.log"
    shell:
        """
        mkdir -p data/raw_external
        wget {params.fasta_path} -O {output.genome_fasta} 
        """


rule get_annotation_gtf:
    message:
        """
        ---Retrieve the annotation in gtf format for a genome.---
        """
    output:
        gtf= "data/raw_external/{genome_name}.gtf.gz"
    params:
        gtf_path = lambda wildcards: config["genomes"][wildcards.genome_name]["gtf"]
    log:
        "logs/{genome_name}_gtf_download.log"
    shell:
        """
        wget {params.gtf_path} -O {output.gtf}
        """


rule get_transcript_fasta:
    message:
        """
        ---Retrieve the transcript sequence in fasta format for a genome.---
        """
    output:
        transcript_fasta="data/raw_external/{genome_name}.transcript.fa.gz"
    params:
        cDNA_path = lambda wildcards: config["genomes"][wildcards.genome_name]["transcript"]
    log:
        "logs/{genome_name}_cdna_download.log"
    shell:
        """
        wget {params.cDNA_path}  -O {output} 
        """
rule index_genome:
    message:
        """
        ---Index a genome using hisat2.---
        """
    input:
        genome_fasta= "data/raw_external/{genome_name}.genome.fa.gz"
    output:
        index_base="hisat_index/{genome_name}",
        index_files_1="hisat_index/{genome_name}.1.bt2",
        index_files_2="hisat_index/{genome_name}.2.bt2",
        index_files_3="hisat_index/{genome_name}.3.bt2",
        index_files_4="hisat_index/{genome_name}.4.bt2",
        index_files_5="hisat_index/{genome_name}.5.bt2",
        index_files_6="hisat_index/{genome_name}.6.bt2",
        index_files_7="hisat_index/{genome_name}.7.bt2",
        index_files_8="hisat_index/{genome_name}.8.bt2",

    params:
        tmp_file=temp("tempfile_{genome_name}.fa")
    log:
        "logs/index__{genome_name}.log"
    shell:
        """
        mkdir -p hisat_index
        # hisat2 cannot use .gz, so unzip to a temporary file first
        gunzip -c {input.genome_fasta} > {params.tmp_file}
        # Build the index
        hisat2-build {params.tmp_file} {output.index_base}/{wildcards.genome_name}
        # Optionally, remove the temporary file (if not automatically handled by Snakemake)
        # rm {params.tmp_file}
        """



rule hisat2_mapping:
    message:
        """
        ---Perform alignment using hisat2.---
        """
    input:
        cleaned_r1="data/cleaned_data/{sample}_1_clean.fq.gz",
        cleaned_r2="data/cleaned_data/{sample}_2_clean.fq.gz",
        index_base=expand("hisat_index/{genome_name}",genome_name=genome_name)
    output:
        bam="aligned_reads/{sample}.bam"
    params:
        extra="--dta"  # example param; include any additional parameters you need
    threads: 6  # adjust the number of threads HISAT2 can use
    log:
        "logs/mapping_{sample}.log"
    shell:
        """
        mkdir -p aligned_reads
        hisat2 -x {input.index_base} \
               -1 {input.trimmed_r1} \
               -2 {input.trimmed_r2} \
               {params.extra} \
               -S /dev/stdout \
               -p {threads} | \
        samtools view -bS - > {output.bam}
        """

rule sort_bam:
    message:
        """
        ---Sorting BAM file.---
        """
    input:
        bam="aligned_reads/{sample}.bam"
    output:
        sorted_bam="sorted_reads/{sample}.sorted.bam"
    log:
        "logs/sorting_{sample}.log"
    threads: 2  # adjust the number of threads samtools can use
    shell:
        """
        mkdir -p sorted_reads
        samtools sort -@ {threads} -o {output.sorted_bam} {input.bam} > {log} 2>&1
        """

rule feature_counts:
    message:
        """
        ---Performing read counting using featureCounts.---
        """
    input:
        sorted_bam=expand("sorted_reads/{sample}.sorted.bam", sample=sample_id),
        gtf= "data/raw_external/{genome_name}.gtf.gz"
    output:
        counts="count_data/{genome_name}/{sample}_counts.txt",
        summary="count_data/{genome_name}/{sample}_counts.txt.summary"
    params:
        extra="-s 2",  # -s 2 specifies that the strandedness is reverse
        tmp_file=temp("tempfile_{genome_name}.gtf")
    threads: 2  # adjust the number of threads featureCounts can use
    log:
        "logs/feature_counts/{genome_name}/{sample}.log"
    shell:
        """
        mkdir -p count_data/{wildcards.genome_name}
        gunzip -c {input.gtf} > {params.tmp_file}
        featureCounts -T {threads} -a {params.tmp_file} -o {output.counts} -p -s 2 {input.sorted_bam} &> {log}
        """


rule salmon_index:
    message:
        """
        ---Creating a Salmon index for the transcriptome.---
        """
    input:
        transcript_fasta="data/raw_external/{genome_name}.transcript.fa.gz"
    output:
        index_dir= directory("salmon_index/{genome_name}")
    params:
        extra="--type quasi"  # Example parameter, adjust as needed
    log:
        "logs/salmon_index/{genome_name}.log"
    threads: 8  # Adjust the number of threads based on your system's capacity
    shell:
        """
        mkdir -p {output.index_dir}
        salmon index -t {input.transcript_fasta} -i {output.index_dir} {params.extra} -p {threads} &> {log}
        """
samples = ["A1", "A2", ...]  # Define a list of sample names

rule salmon_quant:
    message:
        """
        ---Performing transcriptome quantification with Salmon.---
        """
    input:
        index_dir=expand("salmon_index/{genome_name}", genome_name=genome_name),
        cleaned_r1="data/cleaned_data/{sample}_1_clean.fq.gz",
        cleaned_r2="data/cleaned_data/{sample}_2_clean.fq.gz",
    output:
        quant=directory("salmon_quant/{sample}")  # Specify output as a directory
    params:
        extra=""  # Any additional parameters you might want to add
    log:
        "logs/salmon_quant/{sample}.log"
    threads: 6  # adjust the number of threads Salmon can use
    shell:
        """
        mkdir -p {output.quant}
        mkdir -p logs/salmon_quant/{sample}
        salmon quant -i {input.index_dir} \
                     -l ISR \
                     -1 {input.r1} \
                     -2 {input.r2} \
                     -o {output.quant} \
                     -p 6 \
                     {params.extra} \
                     &> {log}
        """
rule remove_rrna:
    input:
        r1="data/trimmed_data/trimmed_{sample}_1.fq.gz",
        r2="data/trimmed_data/trimmed_{sample}_2.fq.gz"
    output:
        rrna_free_r1="data/rrna_removed_data/{sample}_1_rrna_free.fq.gz",
        rrna_free_r2="data/rrna_removed_data/{sample}_2_rrna_free.fq.gz"
    log:
        "logs/remove_rrna_{sample}.log"
    threads: 4
    shell:
        """
        mkdir -p data/rrna_removed_data
        sortmerna --ref rrna/Mus_rrna.fasta,mus_rrna_index --reads1 {input.r1} --reads2 {input.r2} \
                   --aligned rRNA_data/{wildcards.sample}_rRNA \
                   --other data/rrna_removed_data/{wildcards.sample}_rrna_free \
                   -a {threads} &> {log}
        """

rule remove_microbial:
    input:
        rrna_free_r1="data/rrna_removed_data/{sample}_1_rrna_free.fq.gz",
        rrna_free_r2="data/rrna_removed_data/{sample}_2_rrna_free.fq.gz"
    output:
        cleaned_r1="data/cleaned_data/{sample}_1_clean.fq.gz",
        cleaned_r2="data/cleaned_data/{sample}_2_clean.fq.gz"
    log:
        "logs/remove_microbial_{sample}.log"
    threads: 4
    params:
        db_path='kraken2/db'
    shell:
        """
        mkdir -p data/cleaned_data
        # First, classify reads with Kraken2
        kraken2 --db {params.db_path} \
                --paired {input.rrna_free_r1} {input.rrna_free_r2} \
                --threads {threads} \
                --report kraken2_reports/{sample}_report.txt \
                --output kraken2_output/{sample}_output.txt
        
        # Then, filter out unwanted reads based on Kraken2 classification
        # Assume unwanted reads are labeled with 'U' (unclassified) in Kraken2 output
        awk '(NR%4==1){{id=$1}} (NR%4==3){{seq=$1}} (NR%4==0){{qual=$1}} $2!="U"{{print "@"id"\n"seq"\n+\n"qual}}' \
            kraken2_output/{sample}_output.txt > data/cleaned_data/{sample}_clean.fq
        """
