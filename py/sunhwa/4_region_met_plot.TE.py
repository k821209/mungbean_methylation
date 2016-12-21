#!/usr/bin/python3
import sys,pickle
import time
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy.stats import bernoulli, poisson, binom
from scipy.stats import gaussian_kde
from scipy.interpolate import spline



#plt.show()
sns.set_palette("deep", desat=.6)
sns.set_context(rc={"figure.figsize": (8, 4)})

def write_array2file(dic,Outfile_name): 
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()

def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)

data_CG = open_array('S_total.CG.dict.TE_LTR_Copia.gff.100.npy')
data_CHG = open_array('S_total.CHG.dict.TE_LTR_Copia.gff.100.npy')
data_CHH = open_array('S_total.CHH.dict.TE_LTR_Copia.gff.100.npy')
#print(data_CG)
#plt.ylim(0,0.1)
#plt.boxplot(data_CHG)


x = np.arange(150)
print(data_CG)
xnew = np.linspace(x.min(),x.max(),500)
y_CG = np.average(data_CG.T,axis=1)
#y_CG = y_CG/max(y_CG)
sy_CG = spline(x,y_CG,xnew)

y_CHG = np.average(data_CHG.T,axis=1)
#y_CHG = y_CHG/max(y_CHG)
sy_CHG = spline(x,y_CHG,xnew)

y_CHH = np.average(data_CHH.T,axis=1)
#y_CHH = y_CHH/max(y_CHH)
sy_CHH = spline(x,y_CHH,xnew)

plt.plot(xnew,sy_CG,label='CG')
plt.plot(xnew,sy_CHG,label='CHG')
plt.plot(xnew,sy_CHH,label='CHH')

plt.legend()
plt.show()
