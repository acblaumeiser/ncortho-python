import subprocess as sp

def main():
    ncortho = '/media/andreas/Data/ncOrtho/ncortho-python/ncortho/ncortho.py'
    mirnas = '/media/andreas/Data/ncOrtho/sample_data/micrornas/mmu-mir-669a-1.txt'
    models = '/media/andreas/Data/ncOrtho/sample_data/covariance_models'
    output = '/media/andreas/Data/ncOrtho/sample_data/output'
    query = '/media/andreas/Data/ncOrtho/sample_data/genomes/pseudo_genome_for_testing_purposes/pseudo_genome_positive.fa'
    reference = '/media/andreas/Data/ncOrtho/sample_data/genomes/genome/Mus_musculus.chromosomes.fa'
    cpu = 4
    cutoff = '0.6'

    #nc_command = 'python {}'.format(ncortho)
    nc_command = 'python {0} -c {1} -m {2} -n {3} -o {4} -q {5} -r {6} -t {7}'.format(ncortho, cpu, models, mirnas, output, query, reference, cutoff)
    sp.call(nc_command, shell=True)

if __name__ == "__main__":
    main()
