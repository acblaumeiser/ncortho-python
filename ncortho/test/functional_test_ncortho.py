import subprocess as sp

def main():
    cpu = 64
    mirnas = '/home/andreas/Documents/Internship/mouse_project/ortholog_searches/zebrafish/mir-574.tsv'
    models = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_ncortho_search/models'
    ncortho = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/ncortho-python/ncortho/ncortho.py'
    #output = '/home/andreas/Documents/Internship/mouse_project/ortholog_searches/zebrafish/test_anemone'
    output = '/home/andreas/Documents/Internship/mouse_project/ortholog_searches/zebrafish/test'
    #query = '/share/project/andreas/miRNA_project/genomes/zebrafish/Danio_rerio.GRCz11.dna.primary_assembly.fa'
    #query = '/share/project/andreas/miRNA_project/genomes/anemone/Nematostella_vectensis.ASM20922v1.dna.nonchromosomal.fa'
    query = '/home/andreas/Documents/Internship/mouse_project/ortholog_searches/zebrafish/pseudo_genome.fa'
    reference = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/genome/Mus_musculus.chromosomes.fa'
    cutoff = '0.6'
       
    nc_command = 'python {0} -c {1} -m {2} -n {3} -o {4} -q {5} -r {6} -t {7}'.format(ncortho, cpu, models, mirnas, output, query, reference, cutoff)
    #nc_command = 'python {0} -m {2} -n {3} -o {4} -q {5} -r {6} -t {7}'.format(ncortho, cpu, models, mirnas, output, query, reference, cutoff)
    sp.call(nc_command, shell=True)

if __name__ == "__main__":
    main()
