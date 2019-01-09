#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 10:00:22 2018

@author: andreas

Create an miRNA input file for ncOrtho (from miRBase input data)
TODO: include both mature sequences (-5p and -3p) if available
"""

import argparse
from Bio import SeqIO

#convert a FASTA file into a dictionary
def fasta_parser(path):
    with open(path, 'rU') as dictfile:
        seqdict = SeqIO.to_dict(SeqIO.parse(dictfile, 'fasta'))
        return seqdict

#parse a miRBase gtf file and store as a dictionary    
def gtf_parser(gtf, pd, md):
    tmpdict = dict()
    with open(gtf) as gtffile:
        for line in gtffile:
            if not line.startswith('#'):
                                   linedata = line.strip().split()
                                   if linedata[2] == 'miRNA_primary_transcript':
                                       chro = linedata[0]
                                       start = linedata[3]
                                       stop = linedata[4]
                                       strand = linedata[6]
                                       mirnaid = linedata[-1].split('ID=')[1].split(';')[0]
                                       mirnaname = linedata[-1].split('Name=')[1].split(';')[0]
                                       tmpdict[mirnaid] = [mirnaname, chro, start, stop, strand, str(pd[mirnaname].seq), '.', '.']
                                   elif linedata[2] == 'miRNA':
                                       preid = linedata[-1].split('Derives_from=')[1]
                                       matname = linedata[-1].split('Name=')[1].split(';')[0]
                                       try:
                                           matseq = str(md[matname].seq)
                                       except:
                                           print('Sequence for {} not found.'.format(matname))
                                           continue
                                       if '-3p' in matname:
                                           try:
                                               tmpdict[preid][7] = matseq                                  
                                           except:
                                               None
                                       else:
                                           try:
                                               tmpdict[preid][6] = matseq
                                           except:
                                               None                          
    return tmpdict

#write miRNA data into tab separated output file    
def output_writer(out, out_dict):
    with open(out, 'w') as outfile:
        header = '\t'.join(['#miRNA', 'chromosome', 'start', 'end', 'strand', 'pre_seq', 'mat_seq-5p', 'mat_seq-3p']) + '\n'
        outfile.write(header)
        for mirna in out_dict:
                lineout = '\t'.join(out_dict[mirna]) + '\n'
                outfile.write(lineout)
    
def main():
    #Define global variables
    parser = argparse.ArgumentParser(prog='python mirna_input.py', description='tool to produce ncRNA input data files for ncOrtho')
    #GTF file from miRBase
    parser.add_argument('-g', '--gtf', metavar='<.gtf>', type=str, help='path to GTF file')
    #mature sequences
    parser.add_argument('-m', '--mature', metavar='<.fa>', type=str, help='path to mature miRNA sequences')
    #output
    parser.add_argument('-o', '--output', metavar='<path>', type=str, help='path for the output file')
    #precursor sequences
    parser.add_argument('-p', '--pre', metavar='<.fa>', type=str, help='path to pre-miRNA sequences')
    
    args = parser.parse_args()
    
    gtf = args.gtf
    mat = args.mature
    out = args.out
    pre = args.pre
    
    #Pre-computations
    pre_dict = fasta_parser(pre)
    mat_dict = fasta_parser(mat)
    #Output computation
    out_dict = gtf_parser(gtf, pre_dict, mat_dict)
    del pre_dict
    del mat_dict
    #Output writing
    output_writer(out, out_dict)
    
if __name__ == "__main__":
    main()