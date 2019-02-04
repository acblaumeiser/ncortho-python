import unittest
from ... import blastparser

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
        
        #def main():
#    None
#    default_cutoff = 0.8
#    default_msl = 0.9
#    test_mirna = Mirna('2',180389048,180389124,'+','ATCG')
#    samplepath = '/media/andreas/Data/ncOrtho/sample_data/reciproc_blast1.out'
#    
#    bp = BlastParser(test_mirna, samplepath, default_msl)
#    bp.parse_blast_output()
    
#if __name__ == '__main__':
#    main()