#!/usr/bin/env python3
# coding=utf-8

""" Remove non-permissive sites """

"""
Program:    non_permissive
File:       non_permissive.py
Version:    2.0
Date:       09.08.23
Function:   Finds location of non-permissive sequence motifs in specified genome.
Eliminates non-permissive TA sites in genome.

Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/tn_seq

Institution:    Birkbeck University of London
                Project supervisors:  Dr. Irilenia Nobeli
                                  Dr. Sharon Kendall (RVC)
_____________________________________________________________________________
Description:
============
This program uses the sequence motif found to be non-permissive for Himar1 insertions,
(as described by DeJesus et al, 2017) to find non-permissive TA sites in specified
genome and remove the identified sites from the insertion files (.wig).

Usage:
======
non_permissive         SELF

Revision History:
=================
V1.0    24.04.21   Original
V2.0    09.08.23   Updated to iterate through wigs and remove paths

Requirements:
=============
biopython

"""

# **********************************************************************************

def open_ta_list(file):
    """
    Read the .csv file of non-permissive sites and make into list of non-permissive TA sites.

    Input                   file            file of non-perm sites
    Output                  np_sites        list of non-perm sites

    """

    import pandas as pd

    df = pd.read_csv(file, header=1, low_memory=False)
    df = df[['Coordinate', 'Matches Non Permissive Motif']]
    np_sites = []
    for i in range(len(df)):
        if df['Matches Non Permissive Motif'][i] == True:
            np_sites.append(df['Coordinate'][i])

    return (np_sites)


###################################################################################

def remove_np_sites(wig_file, np_sites, outdir):
    """
    Open .wig file of insertion sites, delete non-permissive sites, save new .wig file

    Input               wig_file            file of insertion sites (two columns, position/no insertions)
                        np_sites            python list of positions of non-permissive sites
    Output              new_wig_file        new file of insertion sites (all permissive)

    """
    import pandas as pd
    import os

    # open file and record comment line
    f = open(os.path.expanduser(wig_file), 'r')
    comment = f.readline()
    df = pd.read_csv(wig_file, comment="#", sep=" ", header=0)
    f.close()
    # create name of new file
    new_filename = wig_file.split("/")[-1]
    new_wig_file = outdir + "/perm_" + new_filename
    # open outfile and put in comment
    outfile = open(os.path.expanduser(new_wig_file), 'w')
    outfile.write(comment)
    outfile.write("#non-permissive sites removed" + "\n")
    # find sites not in non-permissable list and write to outfile
    perm_sites = df[~df['variableStep'].isin(np_sites)]
    perm_sites.to_csv(outfile, index=False, sep=' ')
    outfile.close()

    return outfile


###################################################################################

def find_np_sites(motif, offset, fasta):
    """ Use motif finder to find all non-permissive TA sites in specified genome
    Input       motif               motif for non-permissive sites
                offset              number of nt before 'TA' in motif 
                fasta               fasta file for specified genome
    Output      np_sites            list of non-permissive positions

    """
    import os
    from Bio import SeqIO
    from Bio import SeqUtils

    # parse fastq file to extract sequence and convert to string
    fasta_file = os.path.expanduser(fasta)
    for record in SeqIO.parse(fasta_file, "fasta"):
        search_seq = str(record.seq)
    #search_seq = fasta
    # find matches for non-permissive motif
    matches = SeqUtils.nt_search(search_seq.upper(), motif)
    # add one to each position to make index=1 to match positions in wig file
    indexed_positions = []
    for i in range(1, len(matches)):
        new_position = matches[i] + 1 + offset
        indexed_positions.append(new_position)
    # create name of new file
    #np_file = output_file
    # open outfile and put in comment
    #outfile = open(os.path.expanduser(np_file), 'w')
    #outfile.write('\n'.join(map(str, indexed_positions)))
    #outfile.close()
    return indexed_positions

# **********************************************************************************

def iterate_remove_np_sites(wig_dir, outdir, fasta, np_motif='SGNTANCS'):
    """
    Iterate through all .wig files in specified directory, remove sites matching non-permissive motif
    Input       wig_dir             directory of .wig files
                outdir              directory for output files
                fasta               fasta file for specified genome 
                np_motif            non-permissive motif    
    Output      new_file            new file of insertion sites (all permissive)
    """
    
    import glob
    #find all np sites in fastq
    no_sites = find_np_sites(np_motif, 3, fasta)
    wigfiles = glob.glob(wig_dir + "/*_insertions.wig")
    for wig in wigfiles:
        new_file = remove_np_sites(wig, no_sites, outdir)
        print(new_file)


# **********************************************************************************


def test_fasta():

    np_motif_med = 'SGNTANCS'
    #np_motif_long = 'GCGNTANCGC'
    bovis_fasta = "ref_seqs/Mbovis_AF2122_97.fasta"
    #tb_fasta = "ref_seqs/Mtb_H37Rv.fasta"
    # create list of non-permissive positions
    no_sites = find_np_sites(np_motif_med, 3, bovis_fasta)
    print(len(no_sites))
    # remove sites from insertion file, write new file
    new_file = remove_np_sites('tests/output/test_5000_insertions.wig', no_sites, "tests/output/")
    print(new_file)


########## main ###################################################################

if __name__ == "__main__":

    test_fasta()
