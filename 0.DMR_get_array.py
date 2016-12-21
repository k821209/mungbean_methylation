#!/usr/bin/python
'''
scaffold_7      6       -       0       1       CHH     1.0     1       1
scaffold_7      14      +       0       1       CHH     1.0     1       1
scaffold_7      20      -       0       6       CHH     1.0     1       1
scaffold_7      21      +       0       1       CHH     1.0     1       1
scaffold_7      25      -       0       6       CHH     1.0     1       1
scaffold_7      32      +       0       4       CHH     1.0     1       1
scaffold_7      40      -       0       10      CHH     1.0     1       1
scaffold_7      43      -       0       12      CHH     1.0     1       1
scaffold_7      48      +       0       5       CHH     1.0     1       1
scaffold_7      65      +       0       10      CHH     1.0     1       1
scaffold_7      88      +       0       12      CHH     1.0     1       1
'''
from __future__ import print_function
import numpy as np
import sys
sys.path.append('/ref/pipelines/')
import kang
import pickle as pk 
from tqdm import tqdm 
file_list   = ['S_total.qvsig.txt','K_total.qvsig.txt']
file_label  = ['Sunhwa','Kyunggi']
file_ref_fa = 'Vradi_ver6.fa'

dic_ref_fa  = kang.Fasta2dic(file_ref_fa)

bin_size    = 500
threshold   = 10


dicChr2array_M = {}
dicChr2array_U = {}

for chromosome in dic_ref_fa:
    seq   = dic_ref_fa[chromosome]
    dicChr2array_M[chromosome] = np.zeros((len(seq),len(file_list)),dtype=int)
    dicChr2array_U[chromosome] = np.zeros((len(seq),len(file_list)),dtype=int)




# Construct matrix of methylated and unmethylated read numbers for samples
for n,efile in enumerate(file_list):
    for line in tqdm(open(efile)):
        cell   = line.strip().split('\t')
        strChr = cell[0]
        intL   = int(cell[1])-1                          # zerobased
        ## Dont report zero based error for non-ref cultivar 
        #if dic_ref_fa[strChr][intL] == 'C' or dic_ref_fa[strChr][intL] == 'G':
        #    pass
        #else:
        #    print(dic_ref_fa[strChr][intL], dic_ref_fa[strChr][intL-1]) 
        #    print('zerobase error')
        #    exit()
        intM   = int(cell[3])                          # Methylated read number
        intU   = int(cell[4])                          # Un-Methylated read number
        dicChr2array_M[strChr][:,n][intL] += intM
        dicChr2array_U[strChr][:,n][intL] += intU

print("saving matrix..")

pk.dump( dicChr2array_M, open( "dicChr2array_M.p", "wb" ) )
pk.dump( dicChr2array_U, open( "dicChr2array_U.p", "wb" ) )

        





