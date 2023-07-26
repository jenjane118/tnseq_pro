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
        #add barcode to header
        new_name = add_barcode(read1_file, read2_file, output_dir)
        print(new_name)

#********************************************************************************************************

def read_samfile(samfile):
    """
    read sam_file and sort lines between header and reads
    """

    header = []
    reads  = []

    with open(samfile, 'r') as f:
        for line in f:
            line = line.strip()
            if line[0] == "@":
                header.append(line)
            else:
                reads.append(line)
    return header, reads

#********************************************************************************************************

def parse_samflag(read):
    """
    Function to determine strand orientation of alignment from read header
    """
    samflag = read.split()[1]
    if samflag=="16":
        strand = "R"
    elif samflag=="0":
        strand = "F"
    else:
        strand = "*"
    return strand

#********************************************************************************************************

def parse_cigar(cigar):
    """
    Parse cigar string from sam file
    regex: \*|([0-9]+[MIDNSHPX=])+
    Input               cigar               cigar string from sam file
    Output              cig_list            list of cigar string elements
    """
    import re
    #find all number/letter combos in cigar string
    cig_list = re.findall(r'[0-9]+[MIDNSHPX=]+', cigar)
    return cig_list

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

def filter_reads_by_strand(reads):
    """
    Function to filter reads by strand
    Input               reads           list of reads
    Output              fwd_strand_reads    list of reads mapped to forward strand
                        rev_strand_reads    list of reads mapped to reverse strand
    """
    fwd_strand_reads = []
    rev_strand_reads = []
    for read in reads:
        if parse_samflag(read) == "F":
            fwd_strand_reads.append(read)
        else:
            rev_strand_reads.append(read)
    return fwd_strand_reads, rev_strand_reads
    
#**********************************************************************************

def find_ta_position(read, tag_pos):
    """
    Find the TA position relative to + strand
    Input           read                        sequencing read (mapped) from sam file
                    tag_pos                     position of tag in read (from find_tags3)
    Output          ta_position                 predicted position of insertion site in read
                                                (end of tag)


    """
    
    strand      = parse_samflag(read)
    left_pos    = read.split()[3]
    if tag_pos != -1:
        if strand == 'F':
            ta_position = int(left_pos) #(this is a few bases off when sequence matches end of transposon tag)
        else: #reverse strand
            ta_position = int(left_pos) + tag_pos
    else:
        ta_position = "-"
    return ta_position

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
    
    # for forward strand, search space restricted to first 20 nt of read plus length of tag
    seq_space = seq[0:20+len(target_tag)]
    #search string for transposon seq with max num mismatches
    match = mmfind1(seq_space, len(seq_space), target_tag, len(target_tag), max)
    
    return match

#**********************************************************************************

def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.

    from https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value/12141511#12141511
    """
    from bisect import bisect_left
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before
    
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
    with open (fastq_file, 'r') as f:
    # header and sequence for R1 reads
        for index, line in enumerate(f):
            if index % 2 == 0:
                head = line.rstrip()
            # second line and every 4  (sequence)  
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
    # write new file with tagged and no-tagged reads
    new_filename = outdir + "//tag_clipped_" + sample + ".fastq"
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

def remove_dups(fwd_bc_list, rev_bc_list):
    """
    Function to remove duplicate barcodes from fwd and rev lists
    Input               fwd_bc_list         list of fwd barcodes
                        rev_bc_list         list of rev barcodes
    Output              template_bc_list    list of unique barcodes from fwd and rev lists
    
    """
    import numpy as np

    fwd_bc_arr = np.unique(np.array(fwd_bc_list), axis=0)
    rev_bc_arr = np.unique(np.array(rev_bc_list), axis=0)
    fwd_bc_list = fwd_bc_arr.tolist()
    rev_bc_list = rev_bc_arr.tolist()
    #template_bc_list = fwd_bc_list + rev_bc_list
    
    return fwd_bc_list, rev_bc_list

#**********************************************************************************

def remove_dup_reads(sam_file):
    """ 
    Function to remove duplicate reads

    Input               sam_file                            .sam file of mapped reads
                        
    Output              unique_list_f                       list of unique forward reads
                        unique_list_r                       list of unique reverse reads
                                                        
    """
    import numpy as np

    counter = 0
    barcode_list_f = []
    barcode_list_r = []

    #read sam_file and sort lines between header and reads too big for memory
    with open(sam_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("@"):
                pass
            else:
                start_pos = line.split()[3]
                counter += 1
                if counter % 100000 == 0:
                    print("Processed ", counter, "reads")
                #find barcode and strand
                read_name  = line.split()[0]
                strand     = parse_samflag(line)
                barcode    = read_name.split("BC:",1)[1]

                #add to lists of starts and barcodes 
                if strand == "F":
                    barcode_list_f.append([start_pos, barcode])
                else:
                    barcode_list_r.append([start_pos, barcode])

    unique_list_f = list(np.unique(np.array(barcode_list_f), axis=0))
    unique_list_r = list(np.unique(np.array(barcode_list_r), axis=0))

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

def assign_counts_to_sites(ta_sites, template_list_fwd, template_list_rev):
    """
    Function to tally ta sites with insertions
    Input               ta_sites                ordered lists of all possible ta sites in Mbovis genome
                        template_list_fwd       list of lists of barcode, strand, and insertion position for each unique template on fwd strand
                        template_list_rev       list of lists of barcode, strand, and insertion position for each unique template on rev strand          
    Output              read_count_dict         dictionary of ta sites with insertions and read counts

    """

    
    read_count_dict = {}
    #no_match        = []

    fwd_insertions = template_list_fwd
    insertions_list = [line[0] for line in fwd_insertions]
    for site in insertions_list:
        site = int(site)
        #check if ins_site is in list of ta sites (or within 3 nt of ta position)
        closest_ta = take_closest(ta_sites, site)
        if site - closest_ta < 4:
            if closest_ta not in read_count_dict:
                read_count_dict[closest_ta] = 1
            else:
                read_count_dict[closest_ta] += 1
        #else:
            #no_match.append(site)

    rev_insertions = template_list_rev
    insertions_list = [line[0] for line in rev_insertions]
    for site in insertions_list:
        site = int(site)
        closest_ta = take_closest(ta_sites, site)
        if closest_ta - site < 4:  #if alignment includes up to 3 bases of transposon seq will be smaller number than actual ta coord
            if closest_ta not in read_count_dict:
                read_count_dict[closest_ta] = 1
            else:
                read_count_dict[closest_ta] += 1
        #else:
            #no_match.append(site)

    return read_count_dict

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
    fasta       = genome_fasta
    #get genome name from fasta file
    genome      = os.path.basename(fasta).split(".")[0]
    fasta_seq   = open_fasta(genome_fasta)
    ta_sites    = find_insertion_sites(fasta_seq)
    # make list of the positions of every insertion from unique templates
    read_positions = filter_mapped_reads_tag3(samfile, tag="ACTTATCAGCCAACCTGTTA", mismatch_max=2)
    #reduce read positions and barcodes to unique templates
    template_positions = remove_dups(read_positions[0], read_positions[1])
    #count number of reads per ta site
    res=assign_counts_to_sites(ta_sites, template_positions[0], template_positions[1])
    #make complete dictionary of all possible ta sites including those with no insertions
    ta_dict = ta_sites_to_dict(ta_sites, res)
    #write wig file
    wig_file = write_wig_from_dict(ta_dict, sample_name, genome)
    
    return wig_file

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