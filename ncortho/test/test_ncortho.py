#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 12:45:13 2018

@author: andreas
"""

#import unittest
import subprocess

#different start parameters depending on whether running on own pc or workstation
#place = 'ak'
place = 'pc'
    
if place == 'ak':
    cpu = '4'
    mirnas = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/micrornas/mirnas.txt'
    models = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/covariance_models'
    ncortho = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/ncortho-python/ncortho/ncortho.py'
    output = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/output'
    query = '/share/project/andreas/miRNA_project/genomes/NEW/alpaca/Vicugna_pacos.vicPac1.dna.toplevel.fa'
    reference = '/home/andreas/Documents/Internship/M.musculus_root/cm_retry/root/genome/Mus_musculus.chromosomes.fa'
        
elif place == 'pc':
    cpu = '4'
    mirnas = '/media/andreas/Data/ncOrtho/sample_data/micrornas/mmu-mir-669a-1.txt'
    #mirnas = '/media/andreas/Data/ncOrtho/sample_data/micrornas/mirnas.txt'
    models = '/media/andreas/Data/ncOrtho/sample_data/covariance_models'
    ncortho = '/media/andreas/Data/ncOrtho/ncortho-python/ncortho/ncortho.py'
    output = '/media/andreas/Data/ncOrtho/sample_data/output'
    query = '/media/andreas/Data/ncOrtho/sample_data/genomes/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa'
    #query = '/media/andreas/Data/ncOrtho/sample_data/genomes/Vicugna_pacos.vicPac1.dna.toplevel.fa'
    #reference = '/media/andreas/Data/ncOrtho/sample_data/genomes/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa'
    reference = '/media/andreas/Data/ncOrtho/sample_data/genomes/genome/Mus_musculus.chromosomes.fa'
        
nc_command = 'python {0} -c {1} -m {2} -n {3} -o {4} -q {5} -r {6}'.format(ncortho, cpu, models, mirnas, output, query, reference)
#print(nc_command)
subprocess.call(nc_command, shell=True)