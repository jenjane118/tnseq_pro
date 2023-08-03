Jennifer J. Stiens
j.j.stiens@gmail.com
03.08.2023


This is a repository for functions to process Tn-seq data that uses primers designed to multiplex samples and identify PCR duplicates.  The primers were designed similarly to protocol used in Mendum et al, 2019, (BMC Genomics). 

The pipeline for processing the tnseq data from .fastq files to .wig files is explained in the notebook "tnseq_process_pipeline.ipynb" and briefly involves the following processes:

1) Download, quality control and trimming of adapters with fastp
2) Retrieval of unique template barcode from the i7 index reads
3) Trimming the transposon tag from the reads
4) Mapping with BWA-mem
5) Assigning reads to TA sites and removing duplicate PCR products
6) Generation of .wig insertion files

All of the scripts are in the scripts/tnseq_pro.py module and in the snakemake/tnseq directory.
