'''
#Re-implementation of ncOrtho in Python
'''

#import
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
        print('You created a new miRNA object.')

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

def main():
    #Define global variables
    #cpu
    #output
    #input
    #genome of interest
    #covariance model
    #pre-mirna
    #mature mirna
    #coordinates
    #input from text file?
    #parse input parameters
    parser = argparse.ArgumentParser()
    #parser.add_argument('-o', metavar='str', nargs='1')
    cm = False
    output = ''
    mirnas = ''
    query = ''
    msl = 1.0
    mpi = 0

    test = Mirna('chr7',1,10)
    test.information()
######## perform covariance model search ###########
    
    cms_output = ''
    cms_command = ''
    subprocess.call(cms_command, shell=True)
    cm_results = cmsearch_parser(cms_output)
    
    accepted_hits = []
    
    for cmr in cm_results:
        blast_results = blast_search(cmr, r, o)
        bp = BlastParser(None,None,None,None,None)
        bp.parse_blast_output()
        if bp.accepted:
            accepted_hits.append('')
            
    write_output(accepted_hits, o)
    
if __name__ == "__main__":
    main()