#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
sns.set_palette("deep", desat=.6)
sns.set_context(rc={"figure.figsize": (8, 7)})
def head(data,num):
    i = 0 
    for each in data:
            print( each)
            if i == num:
                return()
            i += 1
def get_plot(data,ax,c,l):
	ml=data[:,0]
	bin_num = (max(ml)-min(ml)) / 0.01
	ax.hist(ml, bins=bin_num, color=c, alpha=0.3, label=l)
	#ax.plot(ml,xr)
	ax.legend()
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
file_CG = 'S_total.CG.dict.genemethylv.txt'
file_CHG = 'S_total.CHG.dict.genemethylv.txt'
file_CHH = 'S_total.CHH.dict.genemethylv.txt'
data_CG = np.loadtxt(file_CG,skiprows=1,usecols = (3,4))
data_CHG = np.loadtxt(file_CHG,skiprows=1,usecols = (3,4))
data_CHH = np.loadtxt(file_CHH,skiprows=1,usecols = (3,4))
plt.ylim(0,1000)
plt.xlim(0,1)
get_plot(data_CG,ax1,'r','CG')
get_plot(data_CHG,ax2,'b','CHG')
get_plot(data_CHH,ax3,'g','CHH')

#f.set_size_inches(18.5,10.5)
plt.savefig('mungbean_methyl_sup_fig_7.png',dpi=300)
plt.show()
