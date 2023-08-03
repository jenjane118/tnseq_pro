__author__ = "Jennifer J Stiens"
__copyright__ = "Copyright 2023, Jennifer J Stiens"
__email__ = "j.j.stiens@gmail.com"
__license__ = "MIT"


#snakemake snakefile map with bwa-mem (single end reads for tnseq data), 
#sort with samtools, mapping stats with samtools flagstat

configfile: "config.yaml"

rule all:
    input:
        expand("sorted_reads/flagstat_{sample}.txt", sample=config["samples"]),
        expand("sorted_reads/{sample}.bam.bai", sample=config["samples"]),
        expand("sorted_reads/mapped_{sample}.sam", sample=config["samples"])    

#def get_bwa_map_input_fastqs(wildcards):
 #   return config["samples"][wildcards.sample]

rule bwa_map_single:
    input: 
        fa=config["genomeFile"],
        reads="barcoded_reads/tag_trimmed/tag_clipped_{sample}.reads"
    output: temp("mapped_reads/{sample}.bam")
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'", #this is read-group header if needed
    threads: 8
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    log:
        "logs/bwa_mem/{sample}.log"
    shell: 
        "(bwa mem -t {threads} {input} | "
        "samtools view -Sbh - > {output}) 2> {log}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")
    log:
        "logs/samtools/{sample}.log"
    shell:
        #"(samtools sort -T sorted_reads/{wildcards.sample} "
        #"-O bam {input} > {output}) 2> {log}"
        "(samtools sort -o {output} -O bam {input}) 2> {log}"

rule samtools_filter:
    input: 
        "sorted_reads/{sample}.bam",
        "sorted_reads/{sample}.bam.bai"
    output:
        "sorted_reads/mapped_{sample}.sam"
    log:
        "logs/samtools_filter/{sample}.log"
    shell:
        "(samtools view -h -F 4 {input[0]} -o {output[0]}) 2> {log}" 

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index -b {input}"

rule samtools_flagstat:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/flagstat_{sample}.txt"
    shell:
        "samtools flagstat {input} > {output}"

