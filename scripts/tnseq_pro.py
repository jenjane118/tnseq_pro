#!/usr/bin/env python3

# Copyright 2023
# Jennifer J. Stiens
# j.j.stiens@gmail.com
# 
# This script is part of tnseq_pro read processing software   
# 

#usage: tnseq_pro.py

#**********************************************************************************

def add_barcode(file1, file2, outdir):
    """
    A function to extract the barcode from the associated p7 index read to the header of each read in the fastq file

    Input               file1                       fastq file for sequence reads (R1 or R3)
                        file2                       fastq file for i7 index reads (R2)
                        outdir                      path to output directory
    Output              barcoded<sample>.fastq      new fastq file which has barcode added to header   

    """
    import re
    import os.path

    # identify sample name
    filename = str(file1)
    bn_sample = os.path.basename(file1)
    sample = re.sub(".fastq", "", bn_sample)
    seq_list = []
    with open (file1, 'r') as f1, open(file2, 'r') as f2:
    # header and sequence for R1 reads
        for index, (line1, line2) in enumerate(zip(f1, f2)):
            if index % 4 == 0:
                f1_head = line1.rstrip()
                f2_head = line2.rstrip()
            # second line and every 4  (sequence)  
            if index %4 == 1:
                f1_seq = line1.rstrip()
                #create header and barcode for index reads
                #>A01968:63:H77VYDSX5:4:1101:25455:1423 1:N:0:AACGTGAT BC:GGGGGGGG
                f2_barcode = line2.rstrip()
                # add barcode to end of header in read1 file
                new_head = f1_head + "_BC:" + f2_barcode
                #replace whitespace
                new_head = new_head.replace(" ", "_")
                list_entry = new_head + "\n" + f1_seq + "\n"
                seq_list.append(list_entry)
    # write new file with barcodes
    new_filename = outdir + "/barcode_" + sample + ".fastq"
    with open(new_filename, 'w') as outfile:
        outfile.writelines(seq_list)
    outfile.close()
    return new_filename

#********************************************************************************************************

def iterate_add_barcode(fastq_dir, output_dir):
    """
    Function to iterate through fastq files and add barcode to header.
    Input               fastq_dir           path to directory with fastq files
                        output_dir          path to output directory
    Output                  
    """
    import os
    import glob
    trimmed_files = glob.glob(fastq_dir + "/*_R1_001.fastq")
    print(trimmed_files)
    for read1_file in trimmed_files:
        sample_name = os.path.basename(read1_file)
        read2_file  = "fastq/" + sample_name.replace("R1", "R2")
        read2_file = read2_file.replace("trimmed_", "")
        #read2_file = read2_file.replace(".fastq", ".fastq.gz")
        #add barcode to header
        new_name = add_barcode(read1_file, read2_file, output_dir)
        print(new_name)

#********************************************************************************************************

def mmfind1(G,n,H,m,max): # lengths; assume n>m
  """
  A function to find matches of transposon sequence in read 
  (from tpp_tools.py, https://github.com/mad-lab/transit/blob/master/src/pytpp/tpp_tools.py)
    Input               G                           sequence string
                        n                           length of read sequence
                        H                           pattern string
                        m                           length of pattern
                        max                         maximum mismatches
    Output              i, -1                       start position of match, or -1 for no match
  """

  a = G[:n].find(H[:m])
  if a!=-1: return a # shortcut for perfect matches
  for i in range(0,n-m):  # range of 0 to difference between seq length and len pattern
    cnt = 0
    for k in range(m):
      if G[i+k]!=H[k]: cnt += 1
      if cnt>max: break
    if cnt<=max: return i
  return -1

#********************************************************************************************************

def check_seq(str):
    """
    Is the target tag made of acceptable nucleotide characters
    """
    import re
    chars = set('ACTG')
    seq = set(str.upper())
    if seq.issubset(chars):
        return True
    else:
        return False

#**********************************************************************************

def motif_finder(seq, motif):
    """Finds location of motif (substring) in sequence.
    Input       seq         sequence
                motif       desired subsequence
    Output      locations   list of sequence locations of substrings (one-indexed)
    """
    import re
    locations = []
    p = re.compile('(?='+motif+')')       # use lookahead to find all overlapping matches
    matches = p.finditer(seq)
    for match in matches:
        start = match.span()
        locations.append(start[0] + 1)   # add one to change from zero-indexing to sequence position   


    return locations

#**********************************************************************************

def find_insertion_sites(seq):
    """Finds all the locations of TA insertion sites in sequence

    Only using TA positions from + strand for counting number of sites,
    as they are same site on either strand but start is off by one.
    """
    from Bio.Seq import Seq
    
    seq_obj         = Seq(seq)
    rec_seq         = str(Seq.reverse_complement(seq_obj))
    fwd_positions = motif_finder(seq, "TA")
    rev_positions = motif_finder(rec_seq, "TA")
    
    return fwd_positions

#**********************************************************************************

def open_fasta(refseq):
    with open(refseq, 'r') as file:
        seq = ''
        for line in file:
            if line[0] != '>':
                seq += line.strip()
    return seq

#**********************************************************************************

def find_tags_fastq(seq, target_tag, max):
    
    """
    Find transposon tag in sequence of each forward read. Needs to be within first 20 nt of read (start between 5-10)

      Input           read              mapped read to forward strand (+)
                      target_tag        string that matches transposon sequence
                      max               maximum number of mismatches allowed (default=2)
      Output          start             calculation of start of read/insertion point 
                                        from left-most start position of tag or -1 if no match
    """
    
    # for forward strand, search space restricted to first 22 nt of read plus length of tag (fits len of longest primer)
    seq_space = seq[0:22+len(target_tag)]
    #search string for transposon seq with max num mismatches
    match = mmfind1(seq_space, len(seq_space), target_tag, len(target_tag), max)
    
    return match

#**********************************************************************************

def trim_tag_fastq(fastq_file, outdir, tag="ACTTATCAGCCAACCTGTTA", mismatch_max=2):
    """
    Function to trim transposon tag from fastq reads
    Input               fastq_file          fastq file of barcoded reads
                        target_tag          string of transposon tag
                        max                 maximum number of mismatches allowed
    Output              tag_trimmed_fasta   fasta file of trimmed reads (header and sequence only)
    """
    import os.path
    import re
    filename = str(fastq_file)
    bn_sample = os.path.basename(fastq_file)
    sample = re.sub(".fastq", "", bn_sample)
    tagged_list = []
    notag_list = []
    counter = 0
    with open (fastq_file, 'r') as f:
    # header and sequence for R1 reads
        for index, line in enumerate(f):
            if index % 2 == 0:
                head = line.rstrip()
                counter += 1
                if counter % 5000000 == 0:
                    print("Processed ", counter, "reads")
            # second line and every 2  (sequence)  
            if index % 2 == 1:
                seq = line.rstrip()
                #find transposon tag in new sequence
                tag_start = find_tags_fastq(seq, tag, mismatch_max)
                tag_end = tag_start + len(tag)
                if tag_start != -1:
                    new_seq = seq[tag_end:]
                    tagged_list.append(head + "\n")
                    tagged_list.append(new_seq + "\n")
                else:
                    notag_list.append(head + "\n")
                    notag_list.append(seq + "\n")

    print("Total number of reads processed: ", counter)
    print("Total number of reads with tag: ", len(tagged_list)/2)
    print("Total number of reads without tag: ", len(notag_list)/2)
    print("Percent of no-tagged reads: ", (len(notag_list)/2)/counter)
    
    # write new file with tagged and no-tagged reads
    new_filename = outdir + "/tag_clipped_" + sample + ".fastq"
    with open(new_filename, 'w') as outfile:
        outfile.writelines(tagged_list)
    outfile.close()
    bad_filename = outdir + "/no_tag_" + sample + ".fastq"
    with open(bad_filename, 'w') as outfile:
        outfile.writelines(notag_list)
    
    return new_filename

#**********************************************************************************

def iterate_tag_trim(fastq_directory):
    """
    Function to iterate through fastq files in directory and trim transposon tag
    Input               fastq_directory     directory of fastq files
    Output              
    """
    import glob
    fastq_files = glob.glob(fastq_directory + "/*.fastq")
    for file in fastq_files:
        print("File being processed: ", file)
        trim_tag_fastq(file, "barcoded_reads")

#**********************************************************************************

def remove_dup_reads(sam_file):
    """ 
    Function to remove duplicate reads

    Input               sam_file                            .sam file of mapped reads
                        
    Output              unique_list_f                       list of unique forward reads
                        unique_list_r                       list of unique reverse reads
                                                        
    """
    import numpy as np
    import re

    counter = 0
    barcode_list_f = []
    barcode_list_r = []

    with open(sam_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("@"):
                pass
            else:
                left_pos = int(line.split()[3])
                counter += 1
                if counter % 5000000 == 0:
                    print("Processed ", counter, "reads")
                #find barcode and strand
                read_name  = line.split()[0]
                strand     = int(line.split()[1])
                barcode    = read_name.split("BC:",1)[1]
                cigar      = line.split()[5]
                cig_list   = re.findall(r'[0-9]+[MIDNSHPX=]+', cigar)
                seq_length = int(len(line.split()[9])) #length of sequence

                #add to lists of starts and barcodes 
                if strand == 0:
                    barcode_list_f.append([left_pos, barcode])
                else:
                    if cig_list[0][-1]=="S": #if there is soft-clipping at 3' end of reverse read
                        sc_correction = int(cig_list[0][:-1]) #number of soft-clipped bases
                        rev_start = left_pos + seq_length - sc_correction
                    else:
                        rev_start = left_pos + seq_length
                    barcode_list_r.append([rev_start, barcode])

    unique_list_f = (np.unique(np.array(barcode_list_f), axis=0)).tolist()
    unique_list_r = (np.unique(np.array(barcode_list_r), axis=0)).tolist()

    print("Total number of reads processed: ", counter)
    print("Total number of unique forward reads: ", len(unique_list_f))
    print("Total number of unique reverse reads: ", len(unique_list_r))
    print("Total number of duplicate reads: ", counter - len(unique_list_f) - len(unique_list_r))

    return unique_list_f, unique_list_r

#**********************************************************************************

def ta_sites_to_dict(ta_sites, read_count_dict):
    """ 
    Function to create wig file of insertions
    Input               ta_sites                ordered list of all possible ta sites in Mbovis genome
                        read_count_dict         dictionary of ta sites with insertions and read counts
    Output              wig_dict                dict of insertions

    """
    wig_dict = {}
    for site in ta_sites:
        if site in read_count_dict:
            count = read_count_dict[site]
        else:
            count = 0
        wig_dict[site] = count
    return wig_dict

#**********************************************************************************

def write_wig_from_dict(wig_dict, sample_name, genome):
    """
    Function to write wig file from dict of insertions
    Input               wig_dict            dict of insertions
                        outfile_name        str of outfile name
    Output              outfile             wig file of insertions into output dir
    """
   
    #write template reads
    outfile = "output/" + sample_name + "_insertions.wig"
    with open(outfile, 'w') as f:
        f.write("#generated by tnseq_pro from " + sample_name + "\n")
        f.write("variableStep chrom=" + genome + "\n")
        for key, value in wig_dict.items():
            f.write(f"{key} {value}\n")
    return outfile

#**********************************************************************************

def make_ta_dict(genome_fasta):
    fasta_seq   = open_fasta(genome_fasta)
    ta_sites    = find_insertion_sites(fasta_seq)
    ta_dict     = {}
    for site in ta_sites:
        ta_dict[site] = 0
    return ta_dict

#**********************************************************************************

def assign_counts_to_sites(fasta, template_list_fwd, template_list_rev):
    """
    Function to tally ta sites with insertions
    Input               fasta                   fasta file of genome
                        template_list_fwd       list of lists of barcode and insertion position for each unique template on fwd strand
                        template_list_rev       list of lists of barcode and insertion position for each unique template on rev strand          
    Output              read_count_dict         dictionary of every ta sites with unique read counts

    """

    ta_site_dict = make_ta_dict(fasta)
    no_match        = []

    for line in template_list_fwd:
        site = int(line[0]) -2 #ta starts at -2 from clipped end of tag
        if site in ta_site_dict:
            ta_site_dict[site] += 1
        else:
            no_match.append(line)

    for line in template_list_rev:
        site = int(line[0])
        if site in ta_site_dict:
            ta_site_dict[site] += 1
        else:
            no_match.append(line)

    print("number of unique reads assigned to TA sites: ", sum(ta_site_dict.values()))
    print("number of unique reads with no ta site match: ", len(no_match))

    return ta_site_dict #, no_match

#**********************************************************************************

def sam_to_wig(samfile, genome_fasta, sample_name):
    """
    Wrapper function to create wig file from sam file
    Input           samfile         sam file of mapped reads   
                    genome_fasta    path to fasta file of genome        str
                    sample_name     name of sample                      str
    Output          wig file        wig file of insertions
    """
    import os

    #get genome name from fasta file
    genome      = os.path.basename(genome_fasta).split(".")[0]  
    
    #remove duplicate reads
    unique_reads = remove_dup_reads(samfile)
    fwd_sites = unique_reads[0]
    rev_sites = unique_reads[1]

    #count number of reads per ta site
    res = assign_counts_to_sites(genome_fasta, fwd_sites, rev_sites)
    
    #write wig file
    wig_file = write_wig_from_dict(res, sample_name, genome)
    
    return wig_file

#**********************************************************************************

def iterate_sam_to_wig(sam_dir, genome_fasta):
    # make list of .sam files in directory
    import os
    import glob
    import re
    sam_files = glob.glob(sam_dir + "/*.sam")
    print(sam_files)
    for file in sam_files:
        #find sample name from file
        sample_filename = os.path.basename(file).split(".")[0]
        sample_name = re.findall(r'mapped_(\w*)_R1_001', sample_filename)[0]
        print(sample_name)
        print(file)
        sam_to_wig(file, genome_fasta, sample_name)

#**********************************************************************************

def analyze_dataset(wigfile):
  #can be used to summarise the data in a wig file (from tpp_tools.py)
  data = []
  TAs,ins,reads = 0,0,0
  for line in open(wigfile):
    if line[0]=='#': continue
    if line[:3]=='var': continue # variableStep
    w = line.rstrip().split()
    TAs += 1
    cnt = int(w[1])
    if cnt>1: ins += 1
    reads += cnt
    data.append((cnt,w[0]))

  output = open(wigfile+".stats","w")
  output.write("total TAs: %d, insertions: %d (%0.1f%%), total reads: %d\n" % (TAs,ins,100*(ins/float(TAs)),reads))
  output.write("mean read count per non-zero site: %0.1f\n" % (reads/float(ins)))
  output.write("5 highest counts:\n")
  data.sort(reverse=True)
  for cnt,coord in data[:5]:
    output.write("coord=%s, count=%s\n" % (coord,cnt))
  output.close()