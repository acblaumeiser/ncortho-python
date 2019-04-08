'''
#Re-implementation of ncOrtho in Python
'''

###modules import

###Python
#from __future__ import print_function
import argparse
#import multiprocessing
import os
import subprocess as sp
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

#central class of microRNA objects
class Mirna(object):
    def __init__(self, name, chromosome, start, end, strand, pre, mature, bit):
        self.name = name #miRNA identifier
        self.chromosome = chromosome #chromosome that the miRNA is located on
        self.start = int(start) #start position of the pre-miRNA
        self.end = int(end) #end position of the pre-miRNA
        self.strand = strand #sense (+) or anti-sense (-) strand
        self.pre = pre #nucleotide sequence of the pre-miRNA
        self.mature = mature #nucleotide sequence of the mature miRNA
        #self.star = star #opposite mature sequence, can be functional as well
        #self.mature_5p = mature_5p #nucleotide sequence of the 5p mature miRNA
        #self.mature_3p = mature_3p #nucleotide sequence of the 3p mature miRNA
        self.bit = bit #reference bit score that miRNA receives by its own covariance model

#mirpath: path to file with microRNA data
#cmpath: path to covariance models
#output: outpath for temporary files
def mirna_maker(mirpath, cmpath, output):
    
    mmdict = {} #will be the return value
    
    with open(mirpath) as mirna_file:
        mirna_data = [line.strip().split() for line in mirna_file if not line.startswith('#')]
        #print(mirna_data)
    for mirna in mirna_data:
        mirid = mirna[0]
        # check if the output folder exists, otherwise create it
        if not os.path.isdir('{}/{}'.format(output, mirid)):
            try:
                mkdir = 'mkdir {}/{}'.format(output, mirid)
                sp.call(mkdir, shell=True)
            except:
                print('Cannot create output folder for {}. Skipping to next miRNA.')
                continue
        # obtain the reference bit score for each miRNA by applying it to its own covariance model
        print('Calculating reference bit score for {}.'.format(mirid))
        seq = mirna[5]
        query = '{0}/{1}/{1}.fa'.format(output, mirid)
        model = '{0}/{1}.cm'.format(cmpath, mirid)
        #print(query)
        # check if the covariance model even exists, otherwise skip to the next miRNA
        if not os.path.isfile(model):
            print('No covariance model found for {}.'.format(mirid))
            continue
        
        # create a temporary FASTA file with the miRNA sequence as query for cmsearch
        with open(query, 'w') as tmpfile:
            tmpfile.write('>{0}\n{1}'.format(mirid, seq))
        #cmsearch = 'cmsearch'
        cmsearch = '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch'
        cms_output = '{0}/{1}/cmsearch_{1}_tmp.out'.format(output, mirid)
        cms_log = '{0}/{1}/cmsearch_{1}.log'.format(output, mirid)
        cms_command = '{4} -E 0.01 --noali -o {3} --tblout {0} {1} {2}'.format(cms_output, model, query, cms_log, cmsearch)
        sp.call(cms_command, shell=True)
        with open(cms_output) as cmsfile:
            hits = [line.strip().split() for line in cmsfile if not line.startswith('#')]
            if hits:
                top_score = float(hits[0][14])
            else:
                print('Self bit score not applicable, setting threshold to 0.')
                top_score = 0.0

        mirna.append(top_score)
        rmv_cms = 'rm {}'.format(cms_output)
        rmv_log = 'rm {}'.format(cms_log)
        rmv_fa = 'rm {}'.format(query)
        sp.call(rmv_cms, shell=True)
        sp.call(rmv_log, shell=True)
        sp.call(rmv_fa, shell=True)
        #print(mirna)
        #mmdict[mirna[0]] = Mirna(mirna[0], mirna[1], mirna[2], mirna[3], mirna[4], mirna[5], mirna[6])
        mmdict[mirna[0]] = Mirna(*mirna)
    return mmdict

#cms path to cmsearch output
#cmc cutoff to decide which candidate hits should be included for the reverse BLAST search
#mn name/id of the microRNA
def cmsearch_parser(cms, cmc, mn):
    hits_dict = {}
    chromo_dict = {}
    cut_off = cmc
    print(cut_off)
    #print('!!!!!!!!!!!!!!!!!!!!!!!!')
    with open(cms) as cmsfile:
        hits = [line.strip().split() for line in cmsfile if not line.startswith('#') and float(line.strip().split()[14]) >= cut_off]
        #print(hits)
        if hits:
            #top_score = float(hits[0][14])
            #cut_off = top_score * cmc
            #cut_off = cmc
            #hits_dict = {}
            for i, hit in enumerate(hits):
                #print(hit)
                bit_score = float(hit[14])
                if bit_score >= cut_off:
                    #data = (name, chromosome, start, end, strand)
                    data = ('{0}_c{1}'.format(mn, i+1), hit[0], hit[7], hit[8], hit[9], hit[14])
                    #print(data)
                    #hits_dict[i+1] = data
                    hits_dict[data[0]] = data
#Store the hits that satisfy the bit score cutoff to filter duplicates           
                    try:
                        chromo_dict[data[1]].append(data)
                    except:
                        #None
                        chromo_dict[data[1]] = [data]

                #i += 1
    #print(chromo_dict)
    #return hits_dict
#Loop over the candidate hits to eliminate duplicates

#'ultracontig62': [('mmu-mir-15b_c2', 'ultracontig62', '2169252', '2169306', '+', '49.8'), ('mmu-mir-15b_c3', 'ultracontig62', '2169306', '2169252', '-', '41.5')]
    #print(hits_dict)
    for chromo in chromo_dict:
                nrhits = len(chromo_dict[chromo])
                if nrhits > 1:
                    for hitnr in range(nrhits):
                        start = int(chromo_dict[chromo][hitnr][2])
                        stop = int(chromo_dict[chromo][hitnr][3])
                        strand = chromo_dict[chromo][hitnr][4]
                        score = float(chromo_dict[chromo][hitnr][5])
                        for chitnr in range(hitnr+1, nrhits):
                            if strand != chromo_dict[chromo][chitnr][4]:
                                cstart = int(chromo_dict[chromo][chitnr][2])
                                cstop = int(chromo_dict[chromo][chitnr][3])
                                cscore = float(chromo_dict[chromo][chitnr][5])
    #Test if the two hits from opposite strands overlap, which means one of them is (probably) a false positive
    #Out of two conflicting hits, the one with the highest cmsearch score is retained
                                if start in range(cstart, cstop+1) or stop in range(cstart, cstop+1) or cstart in range(start, stop+1) or cstop in range(start, stop+1):
                                    if score > cscore:
                                        try:
                                            del hits_dict[chromo_dict[chromo][chitnr][0]]
                                        except:
                                            pass
                                            #print(taxon, mirid)
                                            #c += 1
                                    else:
                                        try:
                                            del hits_dict[chromo_dict[chromo][hitnr][0]]
                                        except:
                                            #print(taxon, mirid)
                                            c += 1

    return hits_dict
######################

#['GeneScaffold_587', '-', 'rna_aln', '-', 'cm', '1', '77', '73402', '73326', '-', 'no', '1', '0.31', '0.0', '77.3', '8.3e-16', '!', 'dna:genescaffold', 'genescaffold:vicPac1:GeneScaffold_587:1:293948:1', 'REF']


#s: cmsearch result
#r: reference genome
#o: output name
#c: number of threads
def blast_search(s, r, o, c):
    ### check if BLAST database exists, otherwise create it
    ### database files are .nhr, .nin, .nsq
    ### existence of a file can be checked via os.path.isfile(path)
    file_extensions = ['.nhr', '.nin', '.nsq']
    for fe in file_extensions:
        checkpath = '{}{}'.format(r, fe)
        if not os.path.isfile(checkpath): #at least one of the BLAST db files is not existent and has to be created
            db_command = 'makeblastdb -in {} -dbtype nucl'.format(r)
            print(checkpath)
            print(db_command)
            sp.call(db_command, shell=True)
            break
            
    #blast_output = ''
    blast_command = 'blastn -task blastn -db {0} -query {1} -out {2} -num_threads {3} -outfmt 6'.format(r, s, o, c)
    sp.call(blast_command, shell=True)
    #return blast_output

#a: dictionary of accepted hits
#o: path for output
def write_output(a, o):
    with open(o, 'w') as outfile:
        for hit in a:
            outfile.write('>{0}\n{1}\n'.format(hit, a[hit]))

#def main(ext_args=None):
def main():
    
##### parse command line arguments #####

    #Define global variables
    parser = argparse.ArgumentParser(prog='python ncortho.py', description='ncRNA orthology prediction tool')
    #cpu
    parser.add_argument('-c', '--cpu', metavar='int', type=int, help='number of cpu cores ncOrtho should use')
    #covariance models
    parser.add_argument('-m', '--models', metavar='<path>', type=str, help='path to your covariance models')
    #mirna
    parser.add_argument('-n', '--ncrna', metavar='<path>', type=str, help='path to your reference micrornas')
    #output
    parser.add_argument('-o', '--output', metavar='<path>', type=str, help='path for the output folder')
    #query
    parser.add_argument('-q', '--query', metavar='<.fa>', type=str, help='path to your genome of interest')
    #reference
    parser.add_argument('-r', '--reference', metavar='<.fa>', type=str, help='path to your reference genome')

    args = parser.parse_args()
    #os.getcwd()
    #os.chdir(path)
    #os.path.exists(path)
    #os.path.isfile(path)
    #os.path.isdir(path)

    #cpu = args.cpu        
    ### check if computer provides the desired number of cores
    ### in Python 2 or 3 multiprocessing.cpu_count()
    ### or os.cpu_count() in Python 3

    #cpu = os.cpu_count()
    cpu = args.cpu
    mirnas = args.ncrna
    models = args.models
    output = args.output
    query = args.query
    reference = args.reference
    
    ### default values for testing purposes
    #blast_cutoff = args.blastc
    #cm_cutoff = args.cmc
    #mpi = args.mpi
    #msl = args.msl
    blast_cutoff = 0.8
    cm_cutoff = 0.6
    #mpi = 0
    msl = 0.9
    
    """
    if len(sys.argv) == 1 and not ext_args:
        parser.print_help()
        sys.exit(1)
    elif ext_args:
        args = parser.parse_args(ext_args)
    else:
        args = parser.parse_args()
    """
        
##### create miRNA objects from the list of input miRNAs #####
    
    mirna_dict = mirna_maker(mirnas, models, output)
    #print(mirna_dict)

##### identify candidate orthologs #####
    
    for mir_data in mirna_dict:
        
        ##### perform covariance model search #####

        mirna = mirna_dict[mir_data]
        mirna_id = mirna.name
        outdir = '{}/{}'.format(output, mirna_id)
        #Check if outpath exists, otherwise create it
        if not os.path.isdir(outdir):
            sp.call('mkdir {}'.format(outdir), shell=True)
        #Change to outdir
        os.chdir(outdir)
        print('\n### Running cmsearch for {}. ###\n'.format(mirna_id))
        cms_output = '{0}/cmsearch_{1}.out'.format(outdir, mirna_id)
        cut_off = mirna.bit*cm_cutoff
        #print(cut_off)
        #print(cms_output)
        #infernal = '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch'

        ### original Perl command
        #system("$cmsearch -E 0.01 --cpu $cpu --noali --tblout $cmsearch_out $covariance_model $ukn_genome");
        #cms_command = '{5} -E 0.01 --cpu {0} --noali --tblout {1} {2}/{3}.cm {4}'.format(cpu, cms_output, models, mirna_id, query, infernal)
        cmsearch = '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch'
        #cms_command = '{6} -E 0.01 --incT {5} --cpu {0} --noali --tblout {1} {2}/{3}.cm {4}'.format(cpu, cms_output, models, mirna_id, query, cut_off, cmsearch)
        cms_command = '{6} -T {5} --incT {5} --cpu {0} --noali --tblout {1} {2}/{3}.cm {4}'.format(cpu, cms_output, models, mirna_id, query, cut_off, cmsearch)
        
        #### Remember to include -incT value to limit accepted hits by bit score filter, use reference bit score
        # DONE
        
        #print(cms_command)
        #cms_output = '/media/andreas/Data/ncOrtho/sample_data/output/cmsearch_mmu-mir-1.out'
        #cms_output = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/output/cmsearch_mmu-mir-1.out'

        sp.call(cms_command, shell=True)
        #cm_results = cmsearch_parser(cms_output, cm_cutoff, mirna_id)
        cm_results = cmsearch_parser(cms_output, cut_off, mirna_id)
        #print('!!!!!!!!!!!!!!!')
        print(cm_results)
        
        ##### extract sequences for candidate hits #####
        
        if not cm_results:
            print('\n### No hits found for {}. ###\n'.format(mirna_id))
            continue
    
        else:
            gp = GenomeParser(query, cm_results.values())
        #print(gp.hitlist)
            candidates = gp.extract_sequences()
            #print(candidates)
            #if len(candidates) == 1:
                #print('\n### Search successful, found 1 candidate ortholog. ###\n')
            #else:
            print('\n### Covariance model search successful, found {} ortholog candidate(s). ###\n'.format(len(candidates)))
            
            print('### Evaluating candidates. ###\n')        
        
        ##### perform reverse blast test #####
        
        #accepted_hits = []
        accepted_hits = {}
    
        for candidate in candidates:
            sequence = candidates[candidate]
            temp_fasta = '{0}/{1}.fa'.format(outdir, candidate)
            with open(temp_fasta, 'w') as tempfile:
                tempfile.write('>{0}\n{1}'.format(candidate, sequence))
            blast_output = '{0}/blast_{1}.out'.format(outdir, candidate)
            blast_search(temp_fasta, reference, blast_output, cpu)
            bp = BlastParser(mirna, blast_output, msl)
            #print(bp)
            if bp.parse_blast_output():
                accepted_hits[candidate] = sequence
            #if not bp.parse_blast_output():

                #del candidates[candidate]
                
            #if bp.accepted:
                #accepted_hits.append('')
        #print(accepted_hits)
        ##### write output file #####
        if accepted_hits:
            print('### Writing output of accepted candidates. ###\n')
            outpath = '{0}/{1}_orthologs.fa'.format(outdir, mirna_id)
            write_output(accepted_hits, outpath)
            print('### Finished output writing. ###\n')
        else:
            print('None of the candidates for {} could be verified.\n'.format(mirna_id))
        print('### Finished ortholog search for {}. ###'.format(mirna_id))
if __name__ == "__main__":
    main()
