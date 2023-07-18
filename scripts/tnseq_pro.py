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

#parse samflag 

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

def align_len(read):

    """
    Function to determine read alignment length from start position of reads 
    filtered for transposon tag (should have soft-clipped ~30-33 before gDNA)

    Some cases may exist like 29S89M33S--transposon at both ends? or adapter not trimmed? 
    if forward, can ignore last part of alignment, if reverse, can ignore first part
    
    """
    import re
    cig = read.split()[5]
    cig_list = parse_cigar(cig)
    cig_last = len(cig_list) - 1
    strand = parse_samflag(read)
    if strand == "R":
        # assume alignment without extra bases at 3' end 
        #make sure no soft-clipping at start
        #if 'M' in cig_list[0]:  #this could be 3' end alignment
        if "S" in cig_list[cig_last]:  #check for soft-clipping (tag) at 5' end of read
            match_len = cig_list[cig_last - 1] #aligned part is next part
            #soft_clip = cig_list[cig_last]
            match_len = int(match_len.split("M")[0])
        else:
            match_len = 0
        
    elif strand == "F":
        #for F strand, softclipped bases at 3' end are cigar_list[2] and can be ignored?, if no soft-clipping at 5' end, will be no tag
        match_len = cig_list[1]
        match_len = int(match_len.split("M")[0])
        #soft_clip = cig_list[0]
    else:
        match_len = 0
    return match_len

def rev_start(pos, alignlen):
    """
    Function to find start position for reverse aligned reads
    """
    # need to subtract one from start to match actual start position of read on + strand
    rev_start = int(pos) + alignlen -1
    return rev_start

def parse_cigar(cigar):
    """
    Parse cigar string 
    regex: \*|([0-9]+[MIDNSHPX=])+
    """
    import re
    #keep delimiter
    #cig_list = re.split('([0-9]+M|[0-9]+S)', cigar)
    #remove empty strings
    #cig_list = [i for i in cig_list if i]

    #find all number/letter combos in cigar string
    cig_list = re.findall(r'[0-9]+[MIDNSHPX=]+', cigar)
    return cig_list