import pickle
import time
import numpy as np
from scipy.stats import bernoulli, poisson, binom
import sys

k = np.arange(0, 200)
file_methyl = sys.argv[1]#'/home/k821209/mungbean_methylation/sunhwa/S_1_4_5_1.fastq.gz_bismark_bt2_pe.CX_report.txt.Vr07'
methylC_loc = []
binomial_prob = 0.00111421384126
total_ln = 17096274
i = 0
bj = 0
def write_array2file(dic,Outfile_name): 
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()

def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)
for line in open(file_methyl):
    cell = line.strip().split('\t')
    strS = cell[0]
    #if strS == 'Vr07':
    #    pass
    #else: continue
    strP = cell[1]
    if i == int(total_ln*0.1)*bj:
        print 10*bj, "percent done"
        bj += 1
    i += 1
    strM = cell[3]
    strU = cell[4]
    strT = cell[5]
    fProb_M = round(float(strM) / (float(strM) + float(strU)+0.0001),2)
    site_depth = int(strM) + int(strU)
    rv = binom(site_depth, binomial_prob)
    minimum_number_of_Mreads = k[(rv.pmf(k)<0.05)][0]
    if int(strM) >= minimum_number_of_Mreads:
        #print 'methylated',fProb_M,strM,strU,minimum_number_of_Mreads
        methylC_loc.append(int(strP))
    else:
        #print 'unmethylated',fProb_M,strM,strU,minimum_number_of_Mreads
        pass
methylC_loc = np.array(methylC_loc)
write_array2file(methylC_loc, "/home/k821209/mungbean_methylation/sunhwa/methylC_loc_Vr03.nparray")
