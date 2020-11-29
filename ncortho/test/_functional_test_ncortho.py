import subprocess as sp

def main():
    cpu = 4
    mirnas = '/home/andreas/Documents/Internship/mouse_project/_ortholog_searches/_zebrafish/mir-574.tsv'
    models = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_ncortho_search/models'
    ncortho = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/ncortho-python/ncortho/ncortho.py'
    output = '/home/andreas/Documents/Internship/mouse_project/_ortholog_searches/_zebrafish/test'
    query = '/home/andreas/Documents/Internship/mouse_project/_ortholog_searches/_zebrafish/pseudo_genome.fa'
    reference = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/genome/Mus_musculus.chromosomes.fa'
    cutoff = '0.6'

    #nc_command = 'python {}'.format(ncortho)
    nc_command = 'python {0} -c {1} -m {2} -n {3} -o {4} -q {5} -r {6} -t {7}'.format(ncortho, cpu, models, mirnas, output, query, reference, cutoff)
    sp.call(nc_command, shell=True)

if __name__ == "__main__":
    main()
