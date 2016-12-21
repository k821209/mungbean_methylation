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
sns.set_context(rc={"figure.figsize": (8, 8)})

def write_array2file(dic,Outfile_name): 
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()

def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)

def get_plot(ax,data_CG,data_CHG,data_CHH,t):
	x = np.arange(150)
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
	ax.set_title(t)
#	ax.set_ylim(0,1)
	ax.plot(xnew,sy_CG,label='CG')
	ax.plot(xnew,sy_CHG,label='CHG')
	ax.plot(xnew,sy_CHH,label='CHH')
	ax.legend()
	
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)

data_CG_gpsy  = open_array('S_total.CG.q.dict.TE_LTR_Gypsy.gff.genelimit1k.gff.percent.npy')
data_CHG_gpsy = open_array('S_total.CHG.q.dict.TE_LTR_Gypsy.gff.genelimit1k.gff.percent.npy')
data_CHH_gpsy = open_array('S_total.CHH.q.dict.TE_LTR_Gypsy.gff.genelimit1k.gff.percent.npy')
#print(data_CG)
#plt.ylim(0,0.1)
#plt.boxplot(data_CHG)
data_CG_cpia = open_array('S_total.CG.q.dict.TE_LTR_Copia.gff.genelimit1k.gff.percent.npy')
data_CHG_cpia = open_array('S_total.CHG.q.dict.TE_LTR_Copia.gff.genelimit1k.gff.percent.npy')
data_CHH_cpia = open_array('S_total.CHH.q.dict.TE_LTR_Copia.gff.genelimit1k.gff.percent.npy')

data_CG_tir = open_array('S_total.CG.q.dict.TE_TIR.gff.genelimit1k.gff.percent.npy')
data_CHG_tir = open_array('S_total.CHG.q.dict.TE_TIR.gff.genelimit1k.gff.percent.npy')
data_CHH_tir = open_array('S_total.CHH.q.dict.TE_TIR.gff.genelimit1k.gff.percent.npy')

get_plot(ax1,data_CG_gpsy,data_CHG_gpsy,data_CHH_gpsy,'LTR/Gypsy')
get_plot(ax2,data_CG_cpia,data_CHG_cpia,data_CHH_cpia,'LTR/Copia')
get_plot(ax3,data_CG_tir,data_CHG_tir,data_CHH_tir,'TIR')

plt.savefig('te_methyl.percent.png',dpi=300)
plt.show()
