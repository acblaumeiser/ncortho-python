'''
#Re-implementation of ncOrtho in Python
'''

###modules import

###Python
#from __future__ import print_function
import argparse
import multiprocessing
import os
import subprocess
import sys

###External
#import pyfaidx
#import RNA

###ncOrtho internal modules
from blastparser import BlastParser
from genparser import GenomeParser
#from coreset import CoreSet
#from createcm import CmConstructor

###############################################################################

class Mirna(object):
#central class of microRNA objects
    def __init__(self, name, chromosome, start, end, strand, pre, mature):
        self.name = name #miRNA identifier
        self.chromosome = chromosome #chromosome that the miRNA is located on
        self.start = start #start position of the pre-miRNA
        self.end = end #end position of the pre-miRNA
        self.strand = strand #sense (+) or anti-sense (-) strand
        self.pre = pre #nucleotide sequence of the pre-miRNA
        self.mature = mature #nucleotide sequence of the mature miRNA
        #print('You created a new miRNA object.')

#mirnas path to file with microRNA data
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

#cms path to cmsearch output
#cmc cutoff to decide which candidate hits should be included for the reverse BLAST search
#mn name/id of the microRNA
def cmsearch_parser(cms, cmc, mn):
    hits_dict = {}
    with open(cms) as cmsfile:
        hits = [line.strip().split() for line in cmsfile if not line.startswith('#')]
        if hits:
            top_score = float(hits[0][14])
            cut_off = top_score * cmc
            #hits_dict = {}
            for i, hit in enumerate(hits):
                bit_score = float(hit[14])
                if bit_score >= cut_off:
                    #data = (name, chromosome, start, end, strand)
                    data = ('{0}_c{1}'.format(mn, i+1), hit[0], hit[7], hit[8], hit[9])
                    hits_dict[i+1] = data
        #cmsdict = {}
    return hits_dict

#['GeneScaffold_587', '-', 'rna_aln', '-', 'cm', '1', '77', '73402', '73326', '-', 'no', '1', '0.31', '0.0', '77.3', '8.3e-16', '!', 'dna:genescaffold', 'genescaffold:vicPac1:GeneScaffold_587:1:293948:1', 'REF']


#s: cmsearch result
#r: reference genome
#o: output name
#c: number of threads
def blast_search(s, r, o, c):
    #blast_output = ''
    blast_command = 'blastn -task blastn -db {0} -query {1} -out {2} -num_threads {3} -outfmt 6'.format(r, s, o, c)
    subprocess.call(blast_command, shell=True)
    #return blast_output

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
    #parser = argparse.ArgumentParser()
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
    #output = args.output
    #mirnas = args.mirnas
    #reference = args
    #query = args.query
    #msl = args.msl
    #mpi = args.mpi
    
    ### check if computer provides the desired number of cores
    ### in Python 2 or 3 multiprocessing.cpu_count()
    ### or os.cpu_count() in Python 3
    #cpu = args.cpu
    
    #blast_cutoff = args.blastc
    #cm_cutoff = args.cmc
    
    ### default values for testing purposes
    msl = 0.9
    mpi = 0
    blast_cutoff = 0.8
    cm_cutoff = 0.8
    cpu = os.cpu_count()
    
    #place = 'ak'
    place = 'pc'
    
    if place == 'ak':
        mirnas = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/mirnas.txt'
        models = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/covariance_models'
        output = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/output'
        query = '/share/project/andreas/miRNA_project/genomes/NEW/alpaca/Vicugna_pacos.vicPac1.dna.toplevel.fa'
        reference = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/genome/Mus_musculus.chromosomes.fa'
        
    elif place == 'pc':
        #mirnas = '/media/andreas/Data/ncOrtho/sample_data/micrornas/mmu-mir-669a-1.txt'
        mirnas = '/media/andreas/Data/ncOrtho/sample_data/micrornas/mirnas.txt'
        models = '/media/andreas/Data/ncOrtho/sample_data/covariance_models'
        output = '/media/andreas/Data/ncOrtho/sample_data/output'
        query = '/media/andreas/Data/ncOrtho/sample_data/genomes/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa'
        #query = '/media/andreas/Data/ncOrtho/sample_data/genomes/Vicugna_pacos.vicPac1.dna.toplevel.fa'
        reference = '/media/andreas/Data/ncOrtho/sample_data/genomes/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa'
    
###### create miRNA objects #####
    
    mirna_dict = mirna_maker(mirnas)
    #print(mirna_dict)

##### identify candidate orthologs #####
    
    for mirna in mirna_dict:
        
        ##### perform covariance model search #####

        mirna_id = mirna_dict[mirna].name
        #print(output)
        cms_output = '{0}/cmsearch_{1}.out'.format(output, mirna_id)
        #print(cms_output)
        #infernal = '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch'

        ### original Perl command
        #system("$cmsearch -E 0.01 --cpu $cpu --noali --tblout $cmsearch_out $covariance_model $ukn_genome");
        #cms_command = '{5} -E 0.01 --cpu {0} --noali --tblout {1} {2}/{3}.cm {4}'.format(cpu, cms_output, models, mirna_id, query, infernal)
        cms_command = 'cmsearch -E 0.01 --cpu {0} --noali --tblout {1} {2}/{3}.cm {4}'.format(cpu, cms_output, models, mirna_id, query)
        #print(cms_command)
        #cms_output = '/media/andreas/Data/ncOrtho/sample_data/output/cmsearch_mmu-mir-1.out'
        #cms_output = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/output/cmsearch_mmu-mir-1.out'

        #subprocess.call(cms_command, shell=True)
        cm_results = cmsearch_parser(cms_output, cm_cutoff, mirna_id)
        #print(cm_results)
        
        ##### extract sequences for candidate hits #####
        
        if not cm_results:
            print('No hits for {}.'.format(mirna_id))
            continue
    
        else:
            gp = GenomeParser(query, cm_results.values())
        #print(gp.hitlist)
            candidates = gp.extract_sequences()
            print(candidates)        
        
        ##### perform reverse blast test #####
        
        accepted_hits = []
    
        for candidate in candidates:
            sequence = candidates[candidate]
            temp_fasta = '{0}/{1}.fa'.format(output, candidate)
            with open(temp_fasta, 'w') as tempfile:
                tempfile.write('>{0}\n{1}'.format(candidate, sequence))
            blast_output = '{0}/blast_{1}.out'.format(output, candidate)
            blast_search(temp_fasta, reference, blast_output, cpu)
            bp = BlastParser(mirna_dict[mirna], blast_output, msl)
            bp.parse_blast_output()
            #if bp.accepted:
                #accepted_hits.append('')
        
        ##### write output file #####
        #write_output(accepted_hits, o)
    
if __name__ == "__main__":
    main()