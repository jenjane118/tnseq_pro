__author__ = "Jennifer J Stiens"
__copyright__ = "Copyright 2023, Jennifer J Stiens"
__email__ = "j.j.stiens@gmail.com"
__license__ = "MIT"

#/snakemake/tnseq/fastp/snakefile.smk
#run fastp for filtering and trimming tnseq reads

configfile: "config.yaml"

rule all:
    input:
        expand("trimmed_reads/trimmed_{sample}.fastq.gz", sample=config["samples"])
        
rule fastp:
    input:
        sample="{sample}.fastq.gz"
    output:
        trimmed="trimmed_reads/trimmed_{sample}.fastq.gz",
        failed="trimmed_reads/{sample}.failed_trim.fastq.gz",
        html="fastp/report/{sample}.html",
        json="fastp/report/{sample}.json"
    log:
        "logs/tnseq/fastp/{sample}.log"
    params:
        adapters="-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra=""
    threads: 1
    shell:
        "fastp {params.adapters} -i {input.sample} -o {output.trimmed} "
        "--failed_out {output.failed} -j {output.json} -h {output.html}"
