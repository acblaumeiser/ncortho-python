#Parse the output of a BLAST search and filter the relevant hits

#import Bio
#from __future__ import print_function
#import genparser
#import ncortho

class Mirna(object):
#just for testing purposes, should be inside ncortho main script
    #def __init__(self, chromosome, start, stop, strand, pre, mature):
    def __init__(self, chromosome, start, end, strand, mature):
    #def __init__(self, chromosome, start, end, strand, fivep, threep):
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.mature = mature
        #self.fivep = fivep
        #self.threep = threep
        print('You successfully created a new miRNA object.')

class BlastParser(object):
#object that parses the blast results and tests for the reverse best hit criterion
    def __init__(self, mirna, blastpath, msl):
    #def __init__(self, start, stop, chromosome, strand, blasthits):
        self.start = mirna.start
        self.end = mirna.end
        self.chromosome = mirna.chromosome
        self.strand = mirna.strand
        self.blastpath = blastpath
        self.blasthits = []
        self.msl = msl
        self.top_score = 100
        #self.top_score = blasthits[0][6]
    def __call__(self):
        print('So?')
    def parse_blast_output(self,):
        
        with open(self.blastpath) as blastfile:
            #blasthits = blastfile.readlines()
            blasthits = [line.strip().split() for line in blastfile]
            #print(blasthits)
        for blasthit in blasthits:
            sseqid = blasthit[1]
            sstart = int(blasthit[8])
            send = int(blasthit[9])
            if sseqid == self.chromosome:
                if (sstart <= self.start and self.start <= send) or (sstart <= self.end and self.end <= send):
                    #first within second
                    print("yes")
                    print(blasthit)
                elif (self.start <= sstart and sstart <= self.end) or (self.start <= send and send <= self.end):
                    #second within first
                    print("yeah")
                    print(blasthit)
            ### TODO: calculate length of overlap to see if it is relevant/significant
            ### also the bit score should be compared
            ### add a check for the presence of the mature miRNA

def main():
    None
#    default_cutoff = 0.8
#    default_msl = 0.9
#    test_mirna = Mirna('2',180389048,180389124,'+','ATCG')
#    samplepath = '/media/andreas/Data/ncOrtho/sample_data/reciproc_blast1.out'
#    
#    bp = BlastParser(test_mirna, samplepath, default_msl)
#    bp.parse_blast_output()
    
if __name__ == '__main__':
    main()