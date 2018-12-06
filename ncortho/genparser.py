'''
Extract miRNA sequence from query genome according to cmsearch coordinates
'''

#import Bio
#from Bio import SeqIO
#import pyfaidx

class GenomeParser():
    
    def __init__(self, genome, hits):
        self.genome = genome
        self.hits = hits
    
    def sort_hits(self,):
        #sort hits according to chromosome
        self.hits = sorted(self.hits)
    
    def extract_sequence(self,):
        #identify required chromosomes
        #extract required chromosomes
        required_chromos = set([])
        with open(genome, 'rU') as genomefile:
            for line in genomefile:
                if line.startswith('>'):
                    chromoname = line.strip().split()[0].split('>')[1]
                    if chromoname in required_chromos:
                        sequence = ''
        #dictionary storing the required sequences
        chromosomes = [chromo.seq for chromo in SeqIO.parse(self.genome, 'fasta') if chromo.id in required_chromos] 
        for hit in self.hits:
            #hit of the form (chromosome, start, stop, strand)
            None    

def main():
    genome_parser = GenomeParser()
    genome_parser.sort_hits()
    genome_parser.extract_sequence()
    
if __name__ == "__main__":
    main()