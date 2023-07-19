#!/usr/bin/env python

# Copyright 2023
# Jennifer J. Stiens
# j.j.stiens@gmail.com
# 
# This script is part of tnseq_pro read processing software   

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

def find_ta_position(read, tag_pos):
    """
    Find the TA position relative to + strand
    Input           read                        sequencing read (mapped) from sam file


    """
    strand = parse_samflag(read)
    left_pos = read.split()[3]
    if tag_pos != -1:
        if strand == 'F':
            #TA should be the left most position == soft-clipped bases + tag_pos
            ta_position = int(left_pos)
        else:
            ta_position = int(left_pos) + tag_pos
    else:
        ta_position = "-"
    return ta_position


#**********************************************************************************

def find_tags2(read, target_tag, max):
    """
    Find transposon tag in sequence of each read.

      Input           read              mapped read
                      target_tag        string that matches transposon sequence
                      max               maximum number of mismatches allowed (default=2)
      Output          start             calculation of start of read/insertion point 
                                        from left-most start position of tag or -1 if no match
    """
    from Bio.Seq import Seq
    #find rec of tag
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
    else: #strand is reverse
        match = mmfind1(seq, len(seq), rc_tag, len(rc_tag), max)
        if match != -1:
            start = match  #to get to start of 'TA' in reverse read
        else: 
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