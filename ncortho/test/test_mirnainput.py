# -*- coding: utf-8 -*-

import subprocess

mirna_input = '/media/andreas/Data/ncOrtho/ncortho-python/utils/mirna_input.py'

gtf = '/media/andreas/Data/ncOrtho/sample_data/example/mouse/mmu.gff3'
mat = '/media/andreas/Data/ncOrtho/sample_data/example/mouse/mmu_high_conf_mature.fa'
out = '/media/andreas/Data/ncOrtho/sample_data/example/mouse/mmu_mirna.tsv'
pre = '/media/andreas/Data/ncOrtho/sample_data/example/mouse/mmu_high_conf_hairpin.fa'

command = 'python {0} -g {1} -m {2} -o {3} -p {4}'.format(mirna_input, gtf, mat, out, pre)

subprocess.call(command, shell=True)