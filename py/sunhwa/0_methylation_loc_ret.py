#Vr01    5       -       0       0       CHH     CAT
import pickle
import time
import numpy as np
import scipy
from scipy.stats import bernoulli, poisson, binom
import sys

file_methyl = sys.argv[1]#'/home/k821209/mungbean_methylation/sunhwa/S_1_4_5_1.fastq.gz_bismark_bt2_pe.CX_report.txt.Vr07'
binomial_prob = 0.00111421384126
def write_array2file(dic,Outfile_name): 
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()
def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)
Outfile = open(file_methyl+'.binomed2.txt','w')
for line in open(file_methyl):
    cell = line.strip().split('\t')
    strS = cell[0]
    #if strS == 'Vr07':
    #    pass
    #else: continue
    strP = cell[1]
    strD = cell[2] # strand +/-
    strM = cell[3]
    strU = cell[4]
    strT = cell[5]
    intdp = int(strM) + int(strU) 
    pv = scipy.stats.binom_test(int(strM),intdp,binomial_prob)
    if pv < 0.01:
        print(strS,strP,strD,strM,strU,strT,'1',sep='\t',file=Outfile)
    else:
        print(strS,strP,strD,0,strU,strT,'0',sep='\t',file=Outfile)
Outfile.close()
