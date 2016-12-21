#!/usr/bin/python3

import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy.stats import bernoulli, poisson, binom
from scipy.stats import gaussian_kde
sns.set_palette("deep", desat=.6)
sns.set_context(rc={"figure.figsize": (8, 4)})
def write_array2file(dic,Outfile_name): # arrayµçictµç
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()

def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)
f, (ax1,ax2,ax3) = plt.subplots(3,sharex=True, sharey=True)
file_CG = 'S_total.CG.wincount.npy'
file_CHG = 'S_total.CHG.wincount.npy'
file_CHH = 'S_total.CHH.wincount.npy'
data_CG = open_array(file_CG)
data_CHG = open_array(file_CHG)
data_CHH = open_array(file_CHH)
#print(data_CG)
ax1.hist(data_CG,100, alpha=0.3,color='r',label='CG')
ax1.legend()
ax2.hist(data_CHG,100, alpha=0.3,color='b',label='CHG')
ax2.legend()
ax3.hist(data_CHH,100, alpha=0.3,color='g',label='CHH')
plt.ylim(0,5e6)
ax3.legend()
plt.savefig('mungbean_methyl_sup_fig_8.png',dpi=300)


