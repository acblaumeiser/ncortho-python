'''
#Re-implementation of ncOrtho in Python
'''

#import
#Python
#from __future__ import print_function
import argparse
import multiprocessing
import os
import subprocess
import sys
#import multiprocessing

#import ncOrtho specific modules
from blastparser import BlastParser
from genparser import GenomeParser
#from coreset import CoreSet
#from createcm import CmConstructor

class Mirna(object):
#central class of microRNA objects
    def __init__(self, name, chromosome, start, end, strand, pre, mature):
        self.name = name #miRNA identifier
        self.chromosome = chromosome #chromosome that the miRNA is located on
        self.start = start #start position of the pre-miRNA
        self.end = end #end position of the pre-miRNA
        self.strand = strand #sense or anti-sense strand
        self.pre = pre #nucleotide sequence of the pre-miRNA
        self.mature = mature #nucleotide sequence of the mature miRNA
        #print('You created a new miRNA object.')

def mirna_maker(mirnas):
    mmdict = {}
    with open(mirnas) as mirna_file:
        mirna_data = [line.strip().split() for line in mirna_file if not line.startswith('#')]
        #print(mirna_data)
    for mirna in mirna_data:
        #print(mirna)
        #mmdict[mirna[0]] = Mirna(mirna[0], mirna[1], mirna[2], mirna[3], mirna[4], mirna[5], mirna[6])
        mmdict[mirna[0]] = Mirna(*mirna)
    return mmdict

def cmsearch_parser(cms):
    with open(cms) as cmsfile:
        hits = [line.strip().split() for line in cmsfile if not line.startswith('#')]
        hits_dict = {}
        for i, hit in enumerate(hits):
            #print(i)
            #data = (name, chromosome, start, end, strand)
            data = ('mir-1_c{0}'.format(i+1), hit[0], hit[7], hit[8], hit[9])
            hits_dict[i+1] = data
        #cmsdict = {}
        return hits_dict

#['GeneScaffold_587', '-', 'rna_aln', '-', 'cm', '1', '77', '73402', '73326', '-', 'no', '1', '0.31', '0.0', '77.3', '8.3e-16', '!', 'dna:genescaffold', 'genescaffold:vicPac1:GeneScaffold_587:1:293948:1', 'REF']


#c: cmsearch result
#r: reference genome
#o: output name
def blast_search(c, r, o):
    blast_output = ''
    blast_command = ''
    subprocess.call(blast_command, shell=True)
    return blast_output

#a: list of accepted hits
#o: path for output
def write_output(a, o):
    with open(o, 'w') as outfile:
        for hit in a:
            modified_hit = hit.split()
            outfile.write(modified_hit)

#def main(ext_args=None):
def main():
    
    ##### parse command line arguments #####
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
    #os.getcwd()
    #os.chdir(path)
    #os.path.exists(path)
    #os.path.isfile(path)
    #os.path.isdir(path)
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
    models = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/covariance_models'
    #output = args.output
    output = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/output'
    #mirnas = args.mirnas
    mirnas = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/mirnas.txt'
    #mirnas = '/media/andreas/Data/ncOrtho/sample_data/micrornas/mirnas.txt'
    #reference = args
    reference = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/genome/Mus_musculus.chromosomes.fa'
    #query = args.query
    query = '/share/project/andreas/miRNA_project/genomes/NEW/alpaca/Vicugna_pacos.vicPac1.dna.toplevel.fa'
    #msl = args.msl
    msl = 1.0
    #mpi = args.mpi
    mpi = 0
    #cpu = args.cpu
    ### check if computer provides the desired number of cores
    ### in Python 2 or 3 multiprocessing.cpu_count()
    ### or os.cpu_count() in Python 3
    cpu = 32

    ###### create miRNA objects #####
    
    mirna_dict = mirna_maker(mirnas)
    #print(mirna_dict)

######## perform covariance model search ###########
    
    for mirna in mirna_dict:
        mirna_id = mirna_dict[mirna].name
        #print(output)
        cms_output = '{0}/cmsearch_{1}.out'.format(output, mirna_id)
        #print(cms_output)
        infernal = '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch'
        cms_command = '{5} -E 0.01 --cpu {0} --noali --tblout {1} {2}/{3}.cm {4}'.format(cpu, cms_output, models, mirna_id, query, infernal)
        #print(cms_command)
        #cms_output = '/media/andreas/Data/ncOrtho/sample_data/output/cmsearch_mmu-mir-1.out'
        cms_output = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/output/cmsearch_mmu-mir-1.out'
#system("$cmsearch -E 0.01 --cpu $cpu --noali --tblout $cmsearch_out $covariance_model $ukn_genome");
    subprocess.call(cms_command, shell=True)
    cm_results = cmsearch_parser(cms_output)
    print(cm_results)
    gp = GenomeParser(query, cm_results.values())
    print(gp.hitlist)
    results = gp.extract_sequences()
    print(results)
    
    

##### perform reverse blast test #####
    
    #accepted_hits = []
    
    #for cmr in cm_results:
        #blast_results = blast_search(cmr, r, o)
        #bp = BlastParser(None,None,None,None,None)
        #bp.parse_blast_output()
        #if bp.accepted:
            #accepted_hits.append('')

##### write output file #####            
    #write_output(accepted_hits, o)
    
if __name__ == "__main__":
    main()
