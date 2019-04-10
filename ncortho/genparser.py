'''
ncOrtho submodule
Extract (miRNA) sequence from query genome according to cmsearch coordinates
'''

#import Bio
#from Bio import SeqIO
import pyfaidx
#from pyfaidx import Fasta

class GenomeParser():
    
    def __init__(self, genpath, hitlist):
        self.genpath = genpath
        self.hitlist = hitlist
        self.gene_dict = self.parse_genome()
    
    def parse_genome(self):
        genome = pyfaidx.Fasta(self.genpath)
        return genome
    
    def extract_sequences(self,):
        seq_dict = {}
        for hit in self.hitlist:
            #print('')
            if hit[4] == '+':
                seq = self.gene_dict[hit[1]][int(hit[2])-1:int(hit[3])].seq
            elif hit[4] == '-':
                seq = self.gene_dict[hit[1]][int(hit[3])-1:int(hit[2])].reverse.complement.seq
            seq_dict[hit[0]] = seq
        return seq_dict
        #identify required chromosomes
        #extract required chromosomes
        #required_chromos = set([])
        #with open(self.genome, 'rU') as genomefile:
         #   for line in genomefile:
          #      if line.startswith('>'):
           #         chromoname = line.strip().split()[0].split('>')[1]
            #        if chromoname in required_chromos:
             #           sequence = ''
        #dictionary storing the required sequences
        #chromosomes = [chromo.seq for chromo in SeqIO.parse(self.genome, 'fasta') if chromo.id in required_chromos] 
        #for hit in self.hits:
            #hit of the form (chromosome, start, stop, strand)
            #None    

def main():
    None
'''
    genome_sample = '/media/andreas/Data/ncOrtho/sample_data/genomes/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa'
    #hitlist = '/media/andreas/Data/ncOrtho/sample_data/output/cmsearch_mmu-mir-1_yeast.out'
    hitlist = []
    hitlist.append(('mir-1_c1', 'II', '73326', '73402', '-'))
    hitlist.append(('mir-1_c2', 'VII', '526072', '526147', '+'))
    #print(hitlist)
    gp = GenomeParser(genome_sample, hitlist)
    hits = gp.extract_sequences()
    print(hits)
    #print(gp.hitlist)
    
    
    II     -         rna_aln              -          cm        1       77    73402    73326      -    no    1 0.31   0.0   77.3   8.3e-16 !   dna:genescaffold genescaffold:vicPac1:GeneScaffold_587:1:293948:1 REF
VII    -         rna_aln              -          cm        1       77   526072   526147      +    no    1 0.46   0.0   63.5   1.2e-11 !   dna:genescaffold genescaffold:vicPac1:GeneScaffold_2748:1:649003:1 REF

'''
    #genome_parser = GenomeParser(genome_sample, hitlist)
    #genome_parser.sort_hits()
    #genome_parser.extract_sequence()
    
if __name__ == '__main__':
    main()
