import subprocess as sp

def main():

    #/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/ncortho-python/ncortho/coreset.py
    #script = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/coreset.py'
    script = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/ncortho-python/ncortho/coreset.py'
    ref_gtf_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/ref_gtf'
#core_gtf_paths = glob.glob('/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/core/*.gtf')
    core_gtf_paths = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/core'
    core_fa_paths = '/share/project/andreas/miRNA_project/mouse_core_genomes'
#mirna_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/mmu_mirna.tsv'
#mirna_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/mmu-mir-466c-1.txt'
#mirna_path = sys.argv[1]
#mirna_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/core_set_mirnas.txt'
#mirna_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/mir-411.tsv'
    #mirna_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/mmu_mirna_hairpin_high_confidence.tsv'
    mirna_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/mirnas/test_mirnas.txt'
#reference = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/reference.gtf'
    ref_gtf_path = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/gtf/Mus_musculus.chromosomes.gtf'
    ref_fa_path = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/genome/Mus_musculus.chromosomes.fa'
#reference = '/media/andreas/Data/ncOrtho/sample_data/core_test/gtf/pseudo_ref_genes.gtf'
#reference = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/pseudo_ref_genes.gtf'
#    oma_paths = glob.glob('/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/oma/*')
    oma_paths = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/oma'
    out_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/test_output'
#out_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/output_mir-411'
    #mip = 3
    test_cmd = 'python {} -r {} -c {} -q {} -p {} -o {} -n {}'.format(script, ref_gtf_path, core_gtf_paths, core_fa_paths, oma_paths, out_path, mirna_path)
    sp.call(test_cmd, shell=True)
    #sp.call('python {}'.format(script), shell=True)


if __name__ == "__main__":
    main()
