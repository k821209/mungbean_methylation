#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
sns.set_palette("deep", desat=.6)
sns.set_context("notebook",rc={"figure.figsize": (6, 6)})
def head(data,num):
    i = 0 
    for each in data:
            print( each)
            if i == num:
                return()
            i += 1
def get_plot(data,ax,c,l):
	ml=data[:,0]
	ex=data[:,1]
	(ar,br)=polyfit(ml,ex,1)
	xr=polyval([ar,br],ml)
	ax.scatter(ml,ex,color=c,alpha=0.3,label=l)
	#ax.plot(ml,xr)
	ax.legend()
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
file_CG = 'K_all.CG.q.dict.genemethylv.txt'
file_CHG = 'K_all.CHG.q.dict.genemethylv.txt'
file_CHH = 'K_all.CHH.q.dict.genemethylv.txt'
data_CG = np.loadtxt(file_CG,skiprows=1,usecols = (3,4))
data_CHG = np.loadtxt(file_CHG,skiprows=1,usecols = (3,4))
data_CHH = np.loadtxt(file_CHH,skiprows=1,usecols = (3,4))

plt.xlim(0,1)
plt.ylim(0,100)
get_plot(data_CG,ax1,'r','CG')
get_plot(data_CHG,ax2,'b','CHG')
get_plot(data_CHH,ax3,'g','CHH')

#f.set_size_inches(18.5,10.5)
plt.savefig('kyunggi_met_exp.dotplot.png',dpi=300)
plt.show()
