#!/usr/bin/python
from __future__ import print_function
import numpy as np
import sys
sys.path.append('/ref/pipelines/')
import kang
import pickle as pk
from tqdm import tqdm
from scipy import stats


file_list   = ['S_total.qvsig.txt','K_total.qvsig.txt']
file_label  = ['Sunhwa','Kyunggi']
threshold   = 5
bin_size    = 200
overlap     = 150


file_ref_fa = 'Vradi_ver6.fa'
dic_ref_fa  = kang.Fasta2dic(file_ref_fa)

print('loading pickle..')
dicChr2array_M = pk.load( open( "dicChr2array_M.p", "rb" ) )
dicChr2array_U = pk.load( open( "dicChr2array_U.p", "rb" ) )
print('pickle loaded')



def get_block(array,depth_cut=0):
    lim_len_block = 20
    block_list    = []
    if len(np.shape(array)) == 1:
        rows = 1
        block = []
        for n,j in enumerate(array):
            if j > depth_cut:
                block.append(n)
            else:
                if len(block) > lim_len_block:
                    block_list.append([block[0],block[-1]])
                    block = []
                else:
                    block = []
        if block != []:
            block_list.append([block[0],block[-1]])
    else: 
        rows, columns = np.shape(array)       
        for i in range(rows):
            earray = array[i]
            block = []
            for n,j in enumerate(earray):
                if j > depth_cut:
                    block.append(n)
                else:
                    if len(block) > lim_len_block:
                        block_list.append([i,block[0],block[-1]])
                        block = []
                    else:
                        block = []
    return block_list

Outfile         = open('test.out','w')
Outfile_Mscores = open('Mscores.tst','w')
dicChr2points_A   = {}
dicChr2points_B   = {}

for chromosome in dic_ref_fa.keys():
    #print(chromosome,'processing')
    try:
        M_matrix = dicChr2array_M[chromosome]
        U_matrix = dicChr2array_U[chromosome]
    except KeyError:
        continue
    chromosome_length = len(dic_ref_fa[chromosome]) 
    #bins    = np.arange(0,len(dic_ref_fa[chromosome]),bin_size)
    bins     = [(bin_size-overlap)*x for x in range(chromosome_length/(bin_size-overlap))]
    for n,filename in enumerate(file_list):
        Mscore_list = []
        for ebin in tqdm(bins):
            left  = ebin
            right = ebin+bin_size
            M_array_bin = M_matrix[:,n][left:right]
            U_array_bin = U_matrix[:,n][left:right]
            Mscore      = sum(M_array_bin)
            Mscore_list.append(Mscore)
            print(filename,chromosome,left,right,Mscore,file=Outfile_Mscores)
        Mscore_array = np.array(Mscore_list)
        block_list   = get_block(Mscore_array,depth_cut=threshold)
        for preleft, preright in block_list:
            left  = preleft*(bin_size-overlap)
            right = preright*(bin_size-overlap)+bin_size
            M_array_bin = M_matrix[:,n][left:right]
#            print(left,right,M_array_bin)
            cor_left  = M_array_bin.nonzero()[0][0] + left
            cor_right = M_array_bin.nonzero()[0][-1] + left
            print(chromosome,cor_left,cor_right,file_label[n],sep='\t',file=Outfile)
            if n == 0:
                try:
                    dicChr2points_A[chromosome].append([cor_left,cor_right])
                except KeyError:
                    dicChr2points_A[chromosome] =  [[cor_left,cor_right]]
            elif n == 1:
                try:
                    dicChr2points_B[chromosome].append([cor_left,cor_right])
                except KeyError:
                    dicChr2points_B[chromosome] =  [[cor_left,cor_right]]
           

print('DMR calc..')
Outfile_dmr = open('dmr.out.txt','w')
print('chromosome','pointA','pointB','ref.count.C','ref.count.G','ML_'+file_label[0],'ML_'+file_label[1],file_label[0]+'.Mc',file_label[1]+'.Mc',file_label[0]+'.Uc',file_label[1]+'.Uc','oddsratio', 'pvalue',sep='\t',file=Outfile_dmr)
for chromosome in dicChr2points_A.keys():
    points = dicChr2points_A[chromosome]
    for n,(pointA,pointB) in tqdm(enumerate(points)):
        M_matrix = dicChr2array_M[chromosome]
        U_matrix = dicChr2array_U[chromosome]
        Mcount_A = sum(M_matrix[:,0][pointA:pointB])
        Ucount_A = sum(U_matrix[:,0][pointA:pointB])
        Mcount_B = sum(M_matrix[:,1][pointA:pointB])
        Ucount_B = sum(U_matrix[:,1][pointA:pointB])
        oddsratio, pvalue = stats.fisher_exact([[Mcount_A, Mcount_B], [Ucount_A, Ucount_B]])
        ref_count_C = dic_ref_fa[chromosome][pointA:pointB].count('C')
        ref_count_G = dic_ref_fa[chromosome][pointA:pointB].count('G')
        ML_A     = float(Mcount_A) / (Mcount_A + Ucount_A + 0.000001) 
        ML_B     = float(Mcount_B) / (Mcount_B + Ucount_B + 0.000001)
        if pvalue < 0.05 and ML_A > ML_B:
            print(chromosome,pointA,pointB,ref_count_C,ref_count_G,ML_A,ML_B,Mcount_A, Mcount_B,Ucount_A,Ucount_B,oddsratio, pvalue,sep='\t',file=Outfile_dmr)
    try:
        points = dicChr2points_B[chromosome]
    except KeyError:
        continue
    for n,(pointA,pointB) in tqdm(enumerate(points)):
        M_matrix = dicChr2array_M[chromosome]
        U_matrix = dicChr2array_U[chromosome]
        Mcount_A = sum(M_matrix[:,0][pointA:pointB])
        Ucount_A = sum(U_matrix[:,0][pointA:pointB])
        Mcount_B = sum(M_matrix[:,1][pointA:pointB])
        Ucount_B = sum(U_matrix[:,1][pointA:pointB])
        oddsratio, pvalue = stats.fisher_exact([[Mcount_A, Mcount_B], [Ucount_A, Ucount_B]])
        ref_count_C = dic_ref_fa[chromosome][pointA:pointB].count('C')
        ref_count_G = dic_ref_fa[chromosome][pointA:pointB].count('G')
        ML_A     = float(Mcount_A) / (Mcount_A + Ucount_A + 0.000001) 
        ML_B     = float(Mcount_B) / (Mcount_B + Ucount_B + 0.000001)
        if pvalue < 0.05 and ML_B > ML_A:
            print(chromosome,pointA,pointB,ref_count_C,ref_count_G,ML_A,ML_B,Mcount_A, Mcount_B,Ucount_A,Ucount_B,oddsratio, pvalue,sep='\t',file=Outfile_dmr)
    
        
        


    
        
        


