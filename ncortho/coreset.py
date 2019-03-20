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
import pyfaidx
import sys
#import cPickle as pickle

ref_gtf_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/ref_gtf'
#core_gtf_paths = glob.glob('/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/core/*.gtf')
core_gtf_paths = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/core'
core_fa_paths = '/share/project/andreas/miRNA_project/mouse_core_genomes'
#mirna_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/mmu_mirna.tsv'
#mirna_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/mmu-mir-466c-1.txt'
#mirna_path = sys.argv[1]
mirna_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/core_set_mirnas.txt'
#reference = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/reference.gtf'
reference_path = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/gtf/Mus_musculus.chromosomes.gtf'
#reference = '/media/andreas/Data/ncOrtho/sample_data/core_test/gtf/pseudo_ref_genes.gtf'
#reference = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/pseudo_ref_genes.gtf'
oma_paths = glob.glob('/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/oma/*')
out_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/output'
c = 0
mip = 3

mirna_dict = {}
neighbor_dict = {}

# Parse a GTF file to store the coordinates for each protein-coding gene in a dictionary
def gtf_parser(species):
    species_name = species.split('/')[-1].split('.')[0]
    #print(species)
    #gen_dict = {}
    chr_dict = {}
    chromo = ''

    #with open(inpath) as infile, open(outpath, 'wb') as outfile:
    with open(species) as infile:
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
                #gen_dict[geneid] = (contig, i)
                chr_dict[geneid] = (contig, i)
                try:
                    chr_dict[contig][i] = (geneid, start, end, strand)
                except:
                    chr_dict[contig] = {i: (geneid, start, end, strand)}
                i += 1
    #if gen_dict:
    #    print gen_dict
    #if chr_dict:
    #    print chr_dict
    #print(species_name)
    #print(len(chr_dict.keys()))
    return chr_dict

#Try to find the ortholog for a given reference gene in a core set species
def ortho_search(r_gene):
    #print(r_gene)
    orthologs = {}
    for core_taxon in ortho_dict.keys():
        try:
            ortholog = ortho_dict[core_taxon][r_gene]
            orthologs[core_taxon] = ortholog
                        #print '{0} is the ortholog for {1} in {2}.'.format(left_ortholog, ref_dict[chromo][len(ref_dict[chromo])][0], core_taxon)
            print('{0} is the ortholog for {1} in {2}.'.format(ortholog, r_gene, core_taxon))
            
        except:
                        #print 'No ortholog for {0} found in {1}.'.format(ref_dict[chromo][len(ref_dict[chromo])][0], core_taxon)
            print('No ortholog found for {0} in {1}.'.format(r_gene, core_taxon))
    return orthologs
    #print(ortho_dict)
    #print(orthologs)

#omapaths = ['/home/andreas/Documents/Internship/M.musculus_root/cm_retry/oma/I.tridecemlineatus']

def fasta_parser(taxon, chromosome, start, end, strand):
    seq = ''
    return seq

#reference_gtf = gtf_parser(reference_path)
#if not reference in core_gtf:
#    core_gtf.append(reference)
#data_dict = gtf_parser(core_gtf)

ortho_dict = {}

#Parse the pairwise orthologs
for oma_path in oma_paths:
    taxon = oma_path.split('/')[-1]
    with open(oma_path) as omafile:
    #for line in omafile:
        #print line.split()

#        orthologs = {ref: core for (ref, core) in [(line.split()[0], line.split()[1]) for line in omafile.readlines() if len(line.split()) == 4]}
        orthologs = {ref: core for (ref, core) in [(line.split()[0], line.split()[1]) for line in omafile.readlines()]}
        ortho_dict[taxon] = orthologs

#Read in the miRNA data
with open(mirna_path) as mirfile:
    mirnas = [line.split() for line in mirfile.readlines() if not line.startswith('#')]
    #print mirnas

### mir-1	chr1	51	100	+	GCTTGGGACACATACTTCTTTATATGCCCATATGAACCTGCTAAGCTATGGAATGTAAAGAAGTATGTATTTCAGGC	TGGAATGTAAAGAAGTATGTAT
### chr1	sgd	gene	1	50	.	+	.	gene_id "gene1";

#ref_dict = gtf_parser(reference_path)
#with open(ref_gtf, 'wb') as outfile:
#    pickle.dump(ref_dict, outfile, protocol=2)

with open(ref_gtf_path, 'rb') as tmpfile:
    ref_dict = pickle.load(tmpfile)

#Determine the position of each miRNA and its neighboring gene(s)
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

### case 3): miRNA is located right to the last gene, so the last gene is the left neighbor and there cannot be a right neighbor
    elif start > int(ref_dict[chromo][len(ref_dict[chromo])][2]):
        print('{0} is the left neighbor of {1}.'.format(ref_dict[chromo][len(ref_dict[chromo])][0], mirid))
        print('There is no right neighbor of {0}, since it is located at the end of contig {1}.'.format(mirid, chromo))
        continue

### case 4): miRNA is located either between two genes or overlapping with (an intron of) a gene, either on the same or the opposite strand
    else:
        solved = False
        for i, gene in enumerate(ref_dict[chromo]):
            gene_data = ref_dict[chromo][gene]
            ### case 4.1): miRNA inside gene
            if start >= gene_data[1] and end <= gene_data[2] and strand == gene_data[3]:
                solved = True
                c+=1
                print('{0} is located inside the gene {1}.'.format(mirid, gene_data[0]))
                ortho_hits = ortho_search(gene_data[0])
                for core_tax in ortho_hits:
                    try:
                        neighbor_dict[core_tax][mirid] = ('inside', ortho_hits[core_tax])
                    except:
                        neighbor_dict[core_tax] = {mirid: ('inside', ortho_hits[core_tax])}
                break
            elif start >= gene_data[1] and end <= gene_data[2] and strand != gene_data[3]:
                solved = True
                c+=1
                print('{0} is located opposite of the gene {1}.'.format(mirid, gene_data[0]))
                ortho_hits = ortho_search(gene_data[0])
                for core_tax in ortho_hits:
                    try:
                        neighbor_dict[core_tax][mirid] = ('opposite', ortho_hits[core_tax])
                    except:
                        neighbor_dict[core_tax] = {mirid: ('opposite', ortho_hits[core_tax])}
                break
            elif int(ref_dict[chromo][gene][2]) < start and ref_dict[chromo][gene+1][1] > end:
                solved = True
                #left = gene
                print('{1} is the left neighbor of {2}.'.format(gene, ref_dict[chromo][gene][0], mirid))
                print('{1} is the right neighbor of {2}.'.format(gene, ref_dict[chromo][gene+1][0], mirid))
                #print(gene_data[0])
                #print(ref_dict[chromo][gene+1][0])
                left_hits = ortho_search(gene_data[0])
                right_hits = ortho_search(ref_dict[chromo][gene+1][0])
                #save only the hits where both genes have orthologs in a species
                if left_hits:
                    for taxon in left_hits:
                        if taxon in right_hits:
                            try:
                                neighbor_dict[taxon][mirid] = ('in-between', [left_hits[taxon], right_hits[taxon]])
                            except:
                                neighbor_dict[taxon] = {mirid: ('in-between', [left_hits[taxon], right_hits[taxon]])}
                break                

        if not solved:
            print('Unable to resolve synteny for {}.'.format(mirid))

#print(neighbor_dict)

### Search for the coordinates of the orthologs and extract the sequences
for taxon in neighbor_dict:
    print('Starting synteny analysis for {}'.format(taxon))
    gtf_path = '{0}/{1}.gtf'.format(core_gtf_paths, taxon)
    fasta_path = glob.glob('{0}/{1}*.fa'.format(core_fa_paths, taxon))
    #print(fasta_path)
    #print(fasta_path)
    if len(fasta_path) != 1:
        print('Unable to identify genome file for {}'.format(taxon))
        continue
    genome = pyfaidx.Fasta(fasta_path[0])
    print('Trying to parse GTF file for {}.'.format(taxon))
    try:
        core_gtf_dict = gtf_parser(gtf_path)
        #print(core_gtf_dict)
        #print('Worked')
        #print('Parsed GTF file successfully for {}.'.format(taxon))
        #print(neighbor_dict[taxon])
        for mirna in neighbor_dict[taxon]:
            #print('!!! {} !!!'.format(mirna))
            #print(neighbor_dict[taxon][mirna])
            style = neighbor_dict[taxon][mirna][0]
            #print(style)
            if style == 'inside' or style == 'opposite':
                #print('style')
                #print(style)
                try:
                    ortho_data = core_gtf_dict[neighbor_dict[taxon][mirna][1]]
                    #print(ortho_data)
                    #print(core_gtf_dict[ortho_data[0]][ortho_data[1]])
                    positions = list(core_gtf_dict[ortho_data[0]][ortho_data[1]][1:4])
                    #print(positions)
                    coordinates = [ortho_data[0]] + positions
                    #print(coordinates)
                    seq = genome[coordinates[0]][coordinates[1]-1:coordinates[2]].seq
                    #if coordinates[3] == '+':
                    #    seq = self.gene_dict[hit[1]][int(hit[2])-1:int(hit[3])].seq
                    #elif coordinates[3] == '-':
                    #    seq = self.gene_dict[hit[1]][int(hit[3])-1:int(hit[2])].reverse.complement.seq

                    #print(seq)
                    try:
                        mirna_dict[mirna][taxon] = seq
                    except:
                        mirna_dict[mirna] = {taxon: seq}
                    #print(core_gtf_dict[neighbor_dict[taxon][mirna][1]])
                except:
                    print('{} not found in GTF file.'.format(mirna[1]))
            elif style == 'in-between':
                #print('#\n#\n#\n#\n#\n#\n#\n#\n')
                #print(neighbor_dict[taxon][mirna][2])
                #left_data = core_gtf_dict[neighbor_dict[taxon][mirna][1]]
                left_data = core_gtf_dict[neighbor_dict[taxon][mirna][1][0]]
                #right_data = core_gtf_dict[neighbor_dict[taxon][mirna][2]]
                right_data = core_gtf_dict[neighbor_dict[taxon][mirna][1][1]]
                print('#########################')
#Test to see if the two orthologs are themselves neighbors where their distance cannot be larger than the selected mip value
                if left_data[0] == right_data[0] and abs(left_data[1] - right_data[1]) <= mip:
#Determine which sequence to include for the synteny-based ortholog search depending on the order of orthologs
                    if left_data[1] < right_data[1]:
                        print('left')
                        print(core_gtf_dict[left_data[0]][left_data[1]])
                        print(core_gtf_dict[right_data[0]][right_data[1]])
                        contig = left_data[0]
                        print(contig)
                        seq_start = core_gtf_dict[left_data[0]][left_data[1]][2]
                        print(seq_start)
                        seq_end = core_gtf_dict[right_data[0]][right_data[1]][1]
                        print(seq_end)
                        seq = genome[contig][seq_start-1:seq_end].seq
                        try:
                            mirna_dict[mirna][taxon] = seq
                        except:
                            mirna_dict[mirna] = {taxon: seq}
                    elif right_data[1] < left_data[1]:
                        print('right')
                        print(core_gtf_dict[right_data[0]][right_data[1]])
                        print(core_gtf_dict[left_data[0]][left_data[1]])
                        contig = left_data[0]
                        print(contig)
                        seq_start = core_gtf_dict[left_data[0]][left_data[1]][2]
                        print(seq_start)
                        seq_end = core_gtf_dict[right_data[0]][right_data[1]][1]
                        print(seq_end)
                        seq = genome[contig][seq_start-1:seq_end].seq
                        try:
                            mirna_dict[mirna][taxon] = seq
                        except:
                            mirna_dict[mirna] = {taxon: seq}
                    print('Synteny fulfilled.')
                    #print(left_data)
                    #print(right_data)
                else:
                    print('No shared synteny for {} in {}.'.format(mirna, taxon))
                    print(left_data)
                    print(right_data)
                #print(neighbor_dict[taxon][mirna][1])
                #print(core_gtf_dict[neighbor_dict[taxon][mirna]])
                #print(left_data)
                #print(right_data)
                #print(neighbor_dict[taxon][mirna])
                ### test if the two orthologs are also neighbors
                ### take the end position of the one and the start position of the other to extract sequence
    except:
        #print('No GTF file found for {}'.format(taxon))
        continue
def write_output():
    for mirna in mirna_dict:
        with open('{0}/{1}.fa'.format(out_path, mirna), 'w') as outfile:
            for core_taxon in mirna_dict[mirna]:
                outfile.write('>{0}\n{1}\n'.format(core_taxon, mirna_dict[mirna][core_taxon]))

write_output()

#print(len(mirna_dict))
#print(core_gtf_dict)    
#print(neighbor_dict)
#print(mirna_dict.keys())
#def main():
#    print('coreset')
    
#if __name__ == '__main__':
#    main()
#print(neighbor_dict)
