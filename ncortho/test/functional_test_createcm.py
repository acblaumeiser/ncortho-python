import os
import subprocess as sp

def main():
    createcm_path = (
        '/home/andreas/Documents/Internship/ncOrtho_to_distribute/'
        'ncortho_python/ncortho-python/ncortho/createcm.py'
    )
    test_alignment = (
        '/home/andreas/Documents/Internship/ncOrtho_to_distribute/'
        'ncortho_python/test_covariance_model_construction/mmu-let-7i.sto'
    )
    test_name = test_alignment.split('/')[-1].split('.')[0]
    test_outpath = (
        '/home/andreas/Documents/Internship/ncOrtho_to_distribute/'
        'ncortho_python/test_covariance_model_construction/construction_test/'
        '{}'.format(test_name)
    )
    test_cpu = 4
    create_cmd = 'python {} -a {} -o {} -c {}'.format(
        createcm_path, test_alignment, test_outpath, test_cpu
    )
    sp.call(create_cmd, shell=True)

if __name__ == '__main__':
    main()
