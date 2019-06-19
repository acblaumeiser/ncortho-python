import subprocess as sp

def main():

    #script = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/coreset.py'
    script = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/ncortho-python/ncortho/coreset.py'
    ref_gtf_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/ref_gtf'
    core_gtf_paths = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/core'
    core_fa_paths = '/share/project/andreas/miRNA_project/mouse_core_genomes'
    #mirna_path = '/home/andreas/Documents/Internship/mouse_project/micrornas/mmu_mirna_hairpin_high_confidence.tsv'
    mirna_path = '/home/andreas/Documents/Internship/mouse_project/micrornas/test_set.tsv'
    ref_gtf_path = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/gtf/Mus_musculus.chromosomes.gtf'
    ref_fa_path = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/genome/Mus_musculus.chromosomes.fa'
    oma_paths = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/oma'
    out_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/test_output'
    mgi = 3
    test_cmd = 'python {} -r {} -c {} -q {} -p {} -o {} -n {} -g {} -m {}'.format(script, ref_gtf_path, core_gtf_paths, core_fa_paths, oma_paths, out_path, mirna_path, ref_fa_path, mgi)
    sp.call(test_cmd, shell=True)
    #sp.call('python {}'.format(script), shell=True)


if __name__ == "__main__":
    main()
