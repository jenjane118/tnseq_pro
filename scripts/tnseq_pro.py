#!/usr/bin/env python3

# Copyright 2023
# Jennifer J. Stiens
# j.j.stiens@gmail.com
# 
# This script is part of tnseq_pro read processing software   
# 

#usage: tnseq_pro.py

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
    new_filename = outdir + "barcode_" + sample + ".fastq"
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
    trimmed_files = glob.glob(fastq_dir + "/*.fastq")
    for read1_file in trimmed_files:
        sample_name = os.path.basename(read1_file)
        read2_name = sample_name.replace("_R1_001.fastq", "_R2_001.fastq")
        read2_file = read2_name.replace("trimmed_", "fastq/")
        #add barcode to header
        new_name = add_barcode(read1_file, read2_file, output_dir)
        return new_name

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


#*****************************************************************************************

def filter_mapped_reads3(sam_file, tag="ACTTATCAGCCAACCTGTTA", mismatch_max=2):
  """ 
  Revised filtering function: filters reads for presence of transposon tag and for duplicate barcode/ta-position

  Input             sam_file                            .sam file of mapped reads
                    tag                                 str of transposon tag used, default is Himar1
                    mismatch_max                        integer of maximum mismatches allowed in search for tag
  Output            template barcodes                   list of 3 values: barcode, strand, insert start
                                                        files for template reads, duplicate reads and reads with no transposon tag found
  """
  import sys
  import os
  from operator import itemgetter

  #check tag has valid nucleotides
  if check_seq(tag):
     pass
  else:
     sys.exit("Invalid sequence tag")

  #read sam_file and sort lines between header and reads
  data    = read_samfile(sam_file)
  header  = data[0]
  reads   = data[1]


  barcode_list    = []
  template_reads  = []
  template_barcodes = []
  notags_reads    = []
  dup_reads       = []
  
  for read in reads:
    tag_pos   = find_tags2(read, tag, mismatch_max)
    strand    = parse_samflag(read)
    
    if tag_pos != -1:  #there is match for tranposon tag
        insertion_coord = find_ta_position(read, tag_pos)
        #find barcode
        read_name  = read.split()[0]
        barcode    = read_name.split("BC:",1)[1]
        pos_barcode = [insertion_coord, strand, barcode]
        #add to list of all reads with insertions
        barcode_list.append(pos_barcode) 
        barcode_list.sort(key=itemgetter(0))  #maybe this will speed up search ordering by position--but sorting will add time
        # if hasn't been added before, add read to unique template reads
        if barcode_list.count(pos_barcode) < 2:
          template_reads.append(read)
          template_barcodes.append(pos_barcode)
        else:
          dup_reads.append(read)
    else:
      notags_reads.append(read)
  
  print("Total number of mapped reads: ", len(reads))
  print("number of unique templates (with tag): ", len(template_reads))
  print("number of bad reads (with no tag): ", len(notags_reads))
  print("Number of reads with duplicate barcode/starts: ", len(dup_reads))
  bn = os.path.basename(sam_file)
  #write template reads
  outfile = "tag_filtered_templates_" + bn
  with open(outfile, 'w') as f:
    for line in header:
      f.write(f"{line}\n")
    for line in template_reads:
      f.write(f"{line}\n")
  outfile_dups = "duplicate_reads_" + bn
  with open(outfile_dups, 'w') as f:
    for line in header:
      f.write(f"{line}\n")
    for line in dup_reads:
      f.write(f"{line}\n")
  outfile_notag = "notag_filtered_mapped_" + bn
  with open(outfile_notag, 'w') as f:
    for line in header:
      f.write(f"{line}\n")
    for line in notags_reads:
      f.write(f"{line}\n")

  return template_barcodes

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
    
    return fwd_positions, rev_positions

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

def find_tags_fwd(read, target_tag, max):
    """
    Find transposon tag in sequence of each forward read.

      Input           read              mapped read to forward strand (+)
                      target_tag        string that matches transposon sequence
                      max               maximum number of mismatches allowed (default=2)
      Output          start             calculation of start of read/insertion point 
                                        from left-most start position of tag or -1 if no match
    """
    
    seq = read.split()[9]
    
    #search string for transposon seq with max num mismatches
    match = mmfind1(seq, len(seq), target_tag, len(target_tag), max)
    if match != -1:
        start = match + len(target_tag) #this gives start position of read (start of ta site)
    else:
        start = match
    
    return start

#**********************************************************************************

def find_tags_rev(read, target_tag, max):

    from Bio.Seq import Seq
    #find rc of tag
    tag_seq = Seq(target_tag)
    rc_tag  = str(tag_seq.reverse_complement())
    seq = read.split()[9]
    ##strand is reverse start searching sequence from other end 
    #search reverse sequence (3' to 5') for rc of tag (starting with TAAC--) start position of tag is end of tag and start of insert/ta site
    #in some cases, soft-clipped bases at 3' end of read will make tag too long
    cig = read.split()[5]
    cig_list = parse_cigar(cig)
    if cig_list[0][-1]=="S": #if there is soft-clipping at 3' end of reverse read
        sc_correction = int(cig_list[0][:-1]) #number of soft-clipped bases
        match = mmfind1(seq[sc_correction:len(seq)], len(seq), rc_tag, len(rc_tag), max) #start search from start of aligned bases
    else:
        match = mmfind1(seq, len(seq), rc_tag, len(rc_tag), max)
    
    return match

#**********************************************************************************

def find_tags3(read, target_tag, max):
    """
    Find transposon tag in sequence of each read.

      Input           read              mapped read
                      target_tag        string that matches transposon sequence
                      max               maximum number of mismatches allowed (default=2)
      Output          start             calculation of start of read/insertion point 
                                        from left-most start position of tag or -1 if no match
    """
    from Bio.Seq import Seq
    #find rc of tag
    tag_seq = Seq(target_tag)
    rc_tag  = str(tag_seq.reverse_complement())
    # find strand of read
    strand = parse_samflag(read)
    
    seq = read.split()[9]
    if strand == "F":
        #search string for transposon seq with max num mismatches
        match = mmfind1(seq, len(seq), target_tag, len(target_tag), max)
        if match != -1:
            start = match + len(target_tag) #this gives start position of read (start of ta site)
        else:
            start = match
    else: ##strand is reverse start searching sequence from other end 
        #search reverse sequence (3' to 5') for rc of tag (starting with TAAC--) start position of tag is end of tag and start of insert/ta site
        #in some cases, soft-clipped bases at 3' end of read will make tag too long
        cig = read.split()[5]
        cig_list = parse_cigar(cig)
        if cig_list[0][-1]=="S": #if there is soft-clipping at 3' end of reverse read
            sc_correction = int(cig_list[0][:-1]) #number of soft-clipped bases
            match = mmfind1(seq[sc_correction:len(seq)], len(seq), rc_tag, len(rc_tag), max)
        else:
            match = mmfind1(seq, len(seq), rc_tag, len(rc_tag), max)
        start = match
    
    return start


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

def filter_mapped_reads_tag3(sam_file, tag="ACTTATCAGCCAACCTGTTA", mismatch_max=2):
    """ 
    Revised filtering function: filters reads for presence of transposon tag and for duplicate barcode/ta-position

    Input               sam_file                            .sam file of mapped reads
                        tag                                 str of transposon tag used, default is Himar1
                        mismatch_max                        integer of maximum mismatches allowed in search for tag
    Output              read barcodes                       list of 3 values: insert_start, strand, barcode from reads with insertions
                                                        
    """
    import sys
    import os

    #check tag has valid nucleotides
    if check_seq(tag):
        pass
    else:
        sys.exit("Invalid sequence tag")

    #read sam_file and sort lines between header and reads
    data    = read_samfile(sam_file)
    header  = data[0]
    reads   = data[1]

    #divide reads by strand
    fwd_strand_reads, rev_strand_reads = filter_reads_by_strand(reads)
    
    tagged_reads    = []
    barcode_list_fwd    = []
    barcode_list_rev    = []
    barcode_list_all    = []
    notag_count = 0

    for read in fwd_strand_reads:
        tag_pos   = find_tags_fwd(read, tag, mismatch_max)
        strand    = "F"
        if tag_pos != -1:  #there is match for tranposon tag
            insertion_coord = find_ta_position(read, tag_pos)
            #find barcode
            read_name  = read.split()[0]
            barcode    = read_name.split("BC:",1)[1]
            pos_barcode = [insertion_coord, strand, barcode]
            #add to list of all reads with insertions
            barcode_list_fwd.append(pos_barcode) 
            #append read to tagged reads
            tagged_reads.append(read)
        else:
            notag_count += 1
  
    print("Total number of mapped reads forward strand: ", len(fwd_strand_reads))
    print("number of forward reads with tag: ", len(barcode_list_fwd))
    print("number of forward reads with no tag: ", notag_count)
  
    
    notag_count = 0

    for read in rev_strand_reads:
        tag_pos   = find_tags_rev(read, tag, mismatch_max)
        strand    = "R"
        if tag_pos != -1:  #there is match for tranposon tag
            insertion_coord = find_ta_position(read, tag_pos)
            #find barcode
            read_name  = read.split()[0]
            barcode    = read_name.split("BC:",1)[1]
            pos_barcode = [insertion_coord, strand, barcode]
            #add to list of all reads with insertions
            barcode_list_rev.append(pos_barcode) 
            #append read to tagged reads
            tagged_reads.append(read)
        else:
            notag_count += 1
  
    print("Total number of mapped reads reverse strand: ", len(rev_strand_reads))
    print("number of forward reads with tag: ", len(barcode_list_rev))
    print("number of forward reads with no tag: ", notag_count)

    #write read files (includes all fwd and reverse reads)
    bn = os.path.basename(sam_file)
    outfile = "tag_filtered_templates_" + bn
    with open(outfile, 'w') as f:
        for line in header:
            f.write(f"{line}\n")
        for line in tagged_reads:
            f.write(f"{line}\n")

    #barcode_list_all = barcode_list_fwd + barcode_list_rev
    return barcode_list_fwd, barcode_list_rev

#**********************************************************************************

#function to go through all possible ta sites and create .wig file of insertions

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

def assign_counts_to_sites(ta_sites, barcode_position_lists):
    """
    Function to tally ta sites with insertions
    Input               ta_sites                ordered list of all possible ta sites in Mbovis genome
                        barcode_position_list   list of lists of barcode, strand, and insertion position for each read             
                        strand                  strand of reads in barcode list (F or R)

    """

    
    read_count_dict = {}
    no_match        = []

    fwd_insertions = barcode_position_lists[0]
    insertions_list = [line[0] for line in fwd_insertions]
    for site in insertions_list:
        #check if ins_site is in list of ta sites (or within 3 nt of ta position)
        closest_ta = take_closest(ta_sites, site)
        if site - closest_ta < 4:
            if closest_ta not in read_count_dict:
                read_count_dict[closest_ta] = 1
            else:
                read_count_dict[closest_ta] += 1
        else:
            no_match.append(site)

    rev_insertions = barcode_position_lists[1]
    insertions_list = [line[0] for line in rev_insertions]
    for site in insertions_list:
        closest_ta = take_closest(ta_sites, site)
        if closest_ta - site < 4:  #if alignment includes up to 3 bases of transposon seq will be smaller number than actual ta coord
            if closest_ta not in read_count_dict:
                read_count_dict[closest_ta] = 1
            else:
                read_count_dict[closest_ta] += 1
        else:
            no_match.append(site)

    return read_count_dict, no_match

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
    ta_sites    = find_insertion_sites(fasta_seq)[0]
    # make list of the positions of every insertion from unique templates
    template_positions = filter_mapped_reads_tag3(samfile, tag="ACTTATCAGCCAACCTGTTA", mismatch_max=2)
    #count number of reads per ta site
    res=assign_counts_to_sites(ta_sites, template_positions)
    #make complete dictionary of all possible ta sites including those with no insertions
    ta_dict = ta_sites_to_dict(ta_sites, res[0])
    #write wig file
    wig_file = write_wig_from_dict(ta_dict, sample_name, genome)
    
    return wig_file