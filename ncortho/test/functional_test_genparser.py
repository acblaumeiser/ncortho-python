import os
#import unittest
import sys
sys.path.append(os.path.dirname(os.path.realpath(
    '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python'
    '/ncortho-python/ncortho/genparser.py'))
)
from genparser import GenomeParser

def main():
    sample_genome = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/genome/Mus_musculus.chromosomes.fa'
    #sample_genome = '/share/project/andreas/miRNA_project/genomes/yeast/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa'
    hitlist = [('mmu-mir-574_c1', '2', '1', '78', '+', '107.4')]
    hitlist.append(('mmu-mir-574_c2', '2', '78', '1', '-', '107.4'))
    #hitlist = []
    #hitlist.append(('mir-1_c1', 'II', '2', '50', '+', '50.5'))
    print(hitlist)
    #hitlist.append(('mir-1_c2', 'VII', '526072', '526147', '+'))
    gp = GenomeParser(sample_genome, hitlist)
    hits = gp.extract_sequences()
    print(hits)

#    II     -         rna_aln              -          cm        1       77    73402    73326      -    no    1 0.31   0.0   77.3   8.3e-16 !   dna:genescaffold genescaffold:vicPac1:GeneScaffold_587:1:293948:1 REF
#VII    -         rna_aln              -          cm        1       77   526072   526147      +    no    1 0.46   0.0   63.5   1.2e-11 !   dna:genescaffold genescaffold:vicPac1:GeneScaffold_2748:1:649003:1 REF
   
if __name__ == '__main__':
    main()
