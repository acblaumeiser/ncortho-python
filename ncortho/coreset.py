#Create a core set of orthologs
#Find the corresponding syntenic regions in reference and core species
#Search for core ortholog


#required:
#reference microRNA data (sequence, coordinates)
#reference taxon: genome, blastdb, gtf file with gene coordinates
#core set taxa: genome, gtf file, pairwise orthologs


import glob
import pickle
import genparser
#import pyfaidx
#import cPickle as pickle

#core = glob.glob('/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/core/*.gtf')
ref_gtf_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/ref_gtf'
#core_gtf = glob.glob('/home/andreas/Documents/Internship/M.musculus_root/cm_retry/core/gtf/*.gtf')
core_gtf_paths = glob.glob('/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/core/*.gtf')
core_fa_paths = glob.glob('/share/project/andreas/miRNA_project/mouse_core_genomes/*.fa')
#mirpath = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/mirnas.txt'
#mirpath = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/mmu_mirna.tsv'
#mirpath = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/mirnas/test_mirnas.txt'
mirna_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/mirnas_test_set.tsv'
#reference = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/reference.gtf'
reference_path = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/gtf/Mus_musculus.chromosomes.gtf'
#reference = '/media/andreas/Data/ncOrtho/sample_data/core_test/gtf/pseudo_ref_genes.gtf'
#reference = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/pseudo_ref_genes.gtf'
oma_paths = glob.glob('/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/oma/*')
c=0

def gtf_parser(corespecies):
    species = corespecies.split('/')[-1].split('.')[0]
    print(species)
    gen_dict = {}
    chr_dict = {}
    chromo = ''

    #with open(inpath) as infile, open(outpath, 'wb') as outfile:
    with open(corespecies) as infile:
        for line in infile:
            if not line.startswith('#') and line.split()[2] == 'gene' and line.split('gene_biotype')[1].split('\"')[1] == 'protein_coding':
                linedata = line.strip().split('\t')
                contig = linedata[0]
                geneid = linedata[-1].split('\"')[1]
                start = int(linedata[3])
                end = int(linedata[4])
                strand = linedata[6]
            #print linedata
            #print geneid
                if contig != chromo:
                    i = 1
                    chromo = contig
                gen_dict[geneid] = (contig, i)
                try:
                    chr_dict[contig][i] = (geneid, start, end, strand)
                except:
                    chr_dict[contig] = {i: (geneid, start, end, strand)}
                i += 1
    #if gen_dict:
    #    print gen_dict
    #if chr_dict:
    #    print chr_dict
    return chr_dict

def ortho_search(r_gene):
    orthologs = {}
    for core_taxon in ortho_dict.keys():
        try:
            ortholog = ortho_dict[core_taxon][r_gene]
            orthologs[core_taxon] = ortholog
                        #print '{0} is the ortholog for {1} in {2}.'.format(left_ortholog, ref_dict[chromo][len(ref_dict[chromo])][0], core_taxon)
            print('{0} is the ortholog for {1} in {2}.'.format(ortholog, ref_dict[chromo][gene][0], core_taxon))
            
        except:
                        #print 'No ortholog for {0} found in {1}.'.format(ref_dict[chromo][len(ref_dict[chromo])][0], core_taxon)
            print('No ortholog found for {0} in {1}.'.format(ref_dict[chromo][gene][0], core_taxon))
    return orthologs
    #print(ortho_dict)
    #print(orthologs)

#omapaths = ['/home/andreas/Documents/Internship/M.musculus_root/cm_retry/oma/I.tridecemlineatus']

def fasta_parser(taxon, chromosome, start, end, strand):
    seq = ''
    return seq

reference_gtf = gtf_parser(reference_path)
#if not reference in core_gtf:
#    core_gtf.append(reference)
#data_dict = gtf_parser(core_gtf)

ortho_dict = {}

for oma_path in oma_paths:
    taxon = oma_path.split('/')[-1]
    with open(oma_path) as omafile:
    #for line in omafile:
        #print line.split()

#        orthologs = {ref: core for (ref, core) in [(line.split()[0], line.split()[1]) for line in omafile.readlines() if len(line.split()) == 4]}
        orthologs = {ref: core for (ref, core) in [(line.split()[0], line.split()[1]) for line in omafile.readlines()]}
        ortho_dict[taxon] = orthologs

#print len(orthologs)
#print ortho_dict.keys()

#d = {key: value for (key, value) in iterable}


with open(mirpath) as mirfile:
    mirnas = [line.split() for line in mirfile.readlines() if not line.startswith('#')]
    #print mirnas

### mir-1	chr1	51	100	+	GCTTGGGACACATACTTCTTTATATGCCCATATGAACCTGCTAAGCTATGGAATGTAAAGAAGTATGTATTTCAGGC	TGGAATGTAAAGAAGTATGTAT
### chr1	sgd	gene	1	50	.	+	.	gene_id "gene1";

ref_dict = gtf_parser(reference)
#with open(ref_gtf, 'wb') as outfile:
#    pickle.dump(ref_dict, outfile, protocol=2)

#with open(ref_gtf, 'rb') as tmpfile:
#    ref_dict = pickle.load(tmpfile)



for mirna in mirnas:
    mirid = mirna[0]
    print('### {0} ###'.format(mirid))
    if 'chr' in mirna[1]:
        chromo = mirna[1].split('chr')[1]
    else:
        chromo = mirna[1]
    start = int(mirna[2])
    end = int(mirna[3])
    strand = mirna[4]

### find left neighbor or check if located inside gene
### chr_dict[contig][i] = (geneid, start, end, strand)

### case 1): there is no protein-coding gene on the same contig as the miRNA, so there can be no neighbors (should only occur in highly fragmented assemblies)
    if not chromo in ref_dict.keys():
        print('There are no protein-coding genes on contig {0}. Synteny around {1} cannot be established.'.format(chromo, mirid))
        continue

### case 2): miRNA is located left of the first gene and hence has no left neighbor, the first gene is therefore by default the right neighbor
    if end < int(ref_dict[chromo][1][1]):
        print('There is no left neighbor of {0}, since it is located at the start of contig {1}.'.format(mirid, chromo))
        print('{0} is the right neighbor of {1}.'.format(ref_dict[chromo][1][0], mirid))
        continue
        #found_left = False

### case 3): miRNA is located right to the last gene, so the last gene is the left neighbor and there cannot be a right neighbor
    elif start > int(ref_dict[chromo][len(ref_dict[chromo])][2]):
        #found_left = True
        print('{0} is the left neighbor of {1}.'.format(ref_dict[chromo][len(ref_dict[chromo])][0], mirid))
        print('There is no right neighbor of {0}, since it is located at the end of contig {1}.'.format(mirid, chromo))
        continue

### case 4): miRNA is located either between two genes or inside (an intron of) a gene
    else:
        solved = False
        for i, gene in enumerate(ref_dict[chromo]):
            gene_data = ref_dict[chromo][gene]
            #print gene_data
            ### case 4.1): miRNA inside gene
            if start >= gene_data[1] and end <= gene_data[2] and strand == gene_data[3]:
                solved = True
                c+=1
                print('{0} is located inside the gene {1}.'.format(mirid, gene_data[0]))
                #print(mirna[:5])
                #print(gene_data)
                #try:
                #sequences = {}    
                ortho_hits = ortho_search(gene_data[0])
                ortho_dict[mirid] = ortho_hits
                #for core_taxon in ortho_dict.keys():
                #    try:
                        #print 'Do you even try?'
                 #       left_ortholog = ortho_dict[core_taxon][ref_dict[chromo][len(ref_dict[chromo])][0]]
                        #print '{0} is the ortholog for {1} in {2}.'.format(left_ortholog, ref_dict[chromo][len(ref_dict[chromo])][0], core_taxon)
                  #      print '{0} is the ortholog for {1} in {2}.'.format(left_ortholog, ref_dict[chromo][gene][0], core_taxon)
                   # except:
                        #print 'No ortholog for {0} found in {1}.'.format(ref_dict[chromo][len(ref_dict[chromo])][0], core_taxon)
                    #    print 'No ortholog for {0} found in {1}.'.format(ref_dict[chromo][gene][0], core_taxon)
                break
            elif start >= gene_data[1] and end <= gene_data[2] and strand != gene_data[3]:
                solved = True
                c+=1
                print('{0} is located opposite of the gene {1}.'.format(mirid, gene_data[0]))
                ortho_hits = ortho_search(gene_data[0])
                ortho_dict[mirid] = ortho_hits
                #print(mirna[:5])
                #print(gene_data)
                break
            elif int(ref_dict[chromo][gene][2]) < start and ref_dict[chromo][gene+1][1] > end:
                solved = True
                #left = gene
                print('{1} is the left neighbor of {2}.'.format(gene, ref_dict[chromo][gene][0], mirid))
                print('{1} is the right neighbor of {2}.'.format(gene, ref_dict[chromo][gene+1][0], mirid))
                break
                #ortho_search(gene_data[0])
                #for core_taxon in ortho_dict.keys():
                #    try:
                        #print 'Do you even try?'
                 #       left_ortholog = ortho_dict[core_taxon][ref_dict[chromo][len(ref_dict[chromo])][0]]
                        #print '{0} is the ortholog for {1} in {2}.'.format(left_ortholog, ref_dict[chromo][len(ref_dict[chromo])][0], core_taxon)
                  #      print '{0} is the ortholog for {1} in {2}.'.format(left_ortholog, ref_dict[chromo][gene][0], core_taxon)
                   # except:
                        #print 'No ortholog for {0} found in {1}.'.format(ref_dict[chromo][len(ref_dict[chromo])][0], core_taxon)
                    #    print 'No ortholog for {0} found in {1}.'.format(ref_dict[chromo][gene][0], core_taxon)
        if not solved:
            print('Unable to resolve synteny for {}.'.format(mirid))

### find right neighbor (easy if left neighbor is known)
#1798, ENSMUSG00000101445 is the left neighbor of mmu-mir-1a-1.
    #if found_left:
     #   try:
      #      right = ref_dict[chromo][left+1][0]
       #     print '{1} is the right neighbor of {2}.'.format(left+1, ref_dict[chromo][left+1][0], mirid)
        #except:
         #   print 'No right neighbor for {0}.'.format(mirid)
 #   else:
#        print ref_dict[chromo]
 #       print chromo
    #    if end < ref_dict[chromo][1][1]:
     #       print '{0} is the right neighbor of {1}.'.format(ref_dict[chromo][1][0], mirid)
            
        #print 'No left neighbor for {0}, how about right?'.format(mirid)

#print ref_dict['19']

print(ortho_dict)

class CoreSet:
    def __init__(self):
        None

#def main():
#    print('coreset')
    
#if __name__ == '__main__':
#    main()
