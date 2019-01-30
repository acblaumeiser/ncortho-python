#Parse the output of a BLAST search and filter the relevant hits

#import Bio
#from __future__ import print_function
#import genparser
#import ncortho

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
        if not blasthits:
            #print('Candidate not accepted because the reverse BLAST search was not successful.')
            return False
        else:
            tophit = blasthits[0]
            sseqid = tophit[1]
            sstart = int(tophit[8])
            send = int(tophit[9])
            #print(self.chromosome)
            if sseqid == self.chromosome:
                if (sstart <= self.start and self.start <= send) or (sstart <= self.end and self.end <= send):
                    #first within second
                    #print("yes")
                    #print(tophit)
                    return True
                elif (self.start <= sstart and sstart <= self.end) or (self.start <= send and send <= self.end):
                    #second within first
                    #print("yeah")
                    #print(tophit)
                    return True
#        for blasthit in blasthits:
#            sseqid = blasthit[1]
#            sstart = int(blasthit[8])
#            send = int(blasthit[9])
#            #print(self.chromosome)
#            if sseqid == self.chromosome:
#                if (sstart <= self.start and self.start <= send) or (sstart <= self.end and self.end <= send):
#                    #first within second
#                    print("yes")
#                    print(blasthit)
#                elif (self.start <= sstart and sstart <= self.end) or (self.start <= send and send <= self.end):
#                    #second within first
#                    print("yeah")
#                    print(blasthit)
            
            ### TODO: calculate length of overlap to see if it is relevant/significant
            ### also the bit score should be compared
            ### add a check for the presence of the mature miRNA