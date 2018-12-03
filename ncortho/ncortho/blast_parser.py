#Parse the output of a BLAST search and filter the relevant hits

#import Bio

class BlastParser:
    def __init__(self, start, stop, chromosome, strand, blasthits):
        self.start = start
        self.stop = stop
        self.chromosome = chromosome
        self.strand = strand
        self.blasthits = blasthits
    def __call__(self):
        print 'Now what?'
    def parse_blast_output(self,):
        output = 'No hits found'
        return output
        #check if blast hit corresponds to coordinates of miRNA
        

#def main():
 #   bp = Blastparser()

#if __name__ == "__main__":
  #  main()
