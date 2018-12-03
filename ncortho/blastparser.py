#Parse the output of a BLAST search and filter the relevant hits

#import Bio

class BlastParser(object):
    def __init__(self, mirna, blasthits):
    #def __init__(self, start, stop, chromosome, strand, blasthits):
        self.start = mirna.start
        self.stop = mirna.stop
        self.chromosome = mirna.chromosome
        self.strand = mirna.strand
        self.blasthits = blasthits
        self.top_score = 100
        #self.top_score = blasthits[0][6]
    def __call__(self):
        print 'So?'
    def parse_blast_output(self,):
        for blasthit in blasthits:
            None
            
        print 'No hits found.'
        #output = 'No hits found'
        #return output
        #check if blast hit corresponds to coordinates of miRNA
        
def run_blast_parser():
    bp = BlastParser(1,1,1,1,1)
    return bp
    
#bp = run_blast_parser().parse_blast_output()
#bp()
def main():
    bp = BlastParser(1,1,1,1,1)
    print bp()

if __name__ == "__main__":
    main()