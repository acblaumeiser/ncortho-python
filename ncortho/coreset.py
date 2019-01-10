#Create a core set of orthologs
#Find the corresponding syntenic regions in reference and core species
#Search for core ortholog


#required:
#reference microRNA data (sequence, coordinates)
#reference taxon: genome, blastdb, gtf file with gene coordinates
#core set taxa: genome, gtf file, pairwise orthologs


import glob

core = glob.glob('/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/core/*.gtf')
#core = glob.glob('/home/andreas/Documents/Internship/M.musculus_root/cm_retry/core/gtf/*.gtf')
mirpath = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/mirnas.txt'
#reference = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/reference.gtf'
reference = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/gtf/Mus_musculus.chromosomes.gtf'
#query = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/query.gtf'

core.append(reference)
#print core
gtf_dict = {}

#mirnas = []

#for corespecies in core:
def gtf_parser(corespecies):
    species = corespecies.split('/')[-1].split('.')[0]
    print species
    gen_dict = {}
    chr_dict = {}

    chromo = ''

    #with open(inpath) as infile, open(outpath, 'wb') as outfile:
    with open(corespecies) as infile:
        for line in infile:
            if not line.startswith('#') and line.split()[2] == 'gene':
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

with open(mirpath) as mirfile:
    mirnas = [line.split() for line in mirfile.readlines() if not line.startswith('#')]
    print mirnas

### mir-1	chr1	51	100	+	GCTTGGGACACATACTTCTTTATATGCCCATATGAACCTGCTAAGCTATGGAATGTAAAGAAGTATGTATTTCAGGC	TGGAATGTAAAGAAGTATGTAT
### chr1	sgd	gene	1	50	.	+	.	gene_id "gene1";

ref_dict = gtf_parser(reference)
print ref_dict['14']

for mirna in mirnas:
    mirid = mirna[0]
    chromo = mirna[1]
    start = int(mirna[2])
    end = int(mirna[3])
### find left neighbor
    if start < int(ref_dict[chromo][1][1]):
        print 'No left neighbor for {0}.'.format(mirid)
        found_left = False
    else:
        for gene in ref_dict[chromo]:
            if int(ref_dict[chromo][gene][2]) < start and ref_dict[chromo][gene+1][1] > end:
                found_left = True
                left = gene
                print '{0}, {1} is the left neighbor of {2}.'.format(gene, ref_dict[chromo][gene][0], mirid)
                #print ref_dict[chromo][1538]
                #print ref_dict[chromo][gene+1]

### find right neighbor (easy if left neighbor is known)
#1798, ENSMUSG00000101445 is the left neighbor of mmu-mir-1a-1.
    if found_left:
        right = ref_dict[chromo][left+1][0]
        print '{0}, {1} is the right neighbor of {2}.'.format(left+1, ref_dict[chromo][left+1][0], mirid)
    else:
        print 'No left neighbor for {0}, how about right?'.format(mirid)

#    print ref_dict[chromo]
    

class CoreSet:
    def __init__(self):
        None
        


def main():
    print('coreset')
    
if __name__ == '__main__':
    main()
