'''
#Re-implementation of ncOrtho in Python
'''

#import
from __future__ import print_function
import argparse
import subprocess
import sys

#import ncOrtho specific modules
#import blast_parser
from blastparser import BlastParser
from coreset import CoreSet
from createcm import CmConstructor
from genparser import GenomeParser

class Mirna(object):
    #def __init__(self, chromosome, start, stop, strand, pre, mature):
    def __init__(self, chromosome, start, stop, strand):
        self.chromosome = chromosome
        self.start = start
        self.stop = stop
        self.strand = strand
        print('You created a new miRNA object.')

def mirna_maker(chromosome, start, stop, strand):
    new_mirna = Mirna(chromosome, start, stop, strand)
    return new_mirna

def cmsearch_parser(cms):
    cmsdict = {}
    return cmsdict

#a: list of accepted hits
#o: path for output
def write_output(a, o):
    return None

#c: cmsearch result
#r: reference genome
#o: output name
def blast_search(c, r, o):
    blast_output = ''
    blast_command = ''
    subprocess.call(blast_command, shell=True)
    return blast_output

#def main(ext_args=None):
def main():
    
    #Define global variables
    #cpu
    #output
    #input
    #genome of interest
    #covariance model
    #type of input ncrna, defaul mirna
    #pre-mirna
    #mature mirna
    #coordinates
    #input from text file?
    #parse input parameters
    parser = argparse.ArgumentParser()
    #parser.add_argument('-o', metavar='str', nargs='1', required=True, default='.')
    #...
    
    """
    if len(sys.argv) == 1 and not ext_args:
        parser.print_help()
        sys.exit(1)
    elif ext_args:
        args = parser.parse_args(ext_args)
    else:
        args = parser.parse_args()
    """
    
    #args = parser.parse_args()
    
    #models = args.models
    models = '/media/andreas/Data/ncOrtho/sample_data/covariance_models'
    #output = args.output
    output = '/media/andreas/Data/ncOrtho/sample_data/output'
    #mirnas = args.mirnas
    mirnas = '/media/andreas/Data/ncOrtho/sample_data/micrornas/mirnas.txt'
    #reference = args
    reference = '/media/andreas/Data/ncOrtho/sample_data/genomes/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa'
    #query = args.query
    query = '/media/andreas/Data/ncOrtho/sample_data/genomes/Vicugna_pacos.vicPac1.dna.toplevel.fa'
    #msl = args.msl
    msl = 1.0
    #mpi = args.mpi
    mpi = 0
    #cpu = args.cpu
    cpu = 4
    
    test = Mirna('chr7',1,10,'+')
    #test.information()
######## perform covariance model search ###########
    
    cms_output = ''
    cms_command = ''
    #subprocess.call(cms_command, shell=True)
    #cm_results = cmsearch_parser(cms_output)
    
    accepted_hits = []
    
    #for cmr in cm_results:
        #blast_results = blast_search(cmr, r, o)
        #bp = BlastParser(None,None,None,None,None)
        #bp.parse_blast_output()
        #if bp.accepted:
            #accepted_hits.append('')
            
    #write_output(accepted_hits, o)
    
if __name__ == "__main__":
    main()