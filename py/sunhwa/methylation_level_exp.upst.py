#!/usr/bin/python3 
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

def get_plot(file_in,ax_in,maxml):
	data = np.loadtxt(file_in,usecols=(3,4,),skiprows=1)
	exp = []
	lab = []
	for s in range(maxml):
		bMet = (s*10 <= data[:,0]*100) & (data[:,0]*100 < (s+1)*10)
	#	print(data[:,0])
	#	print(bMet)
		exp_array = data[:,1][bMet]
		exp_array = exp_array[(exp_array<100)]
		exp_array = exp_array
		exp_array = list(exp_array)
		if exp_array == []:
			exp.append([-1])
		else:
			exp.append(exp_array)
		lab.append("%d~%d%%"%(s*10,(s+1)*10))
	lab = np.array(lab)
	exp = np.array(exp)
	sns.boxplot(exp,names=lab, ax=ax_in,widths=0.5)
	
f, (ax1, ax2, ax3) = plt.subplots(3,sharex=True)
get_plot('S_total.CG.dict.upstream.methylv.txt',ax1,10)
get_plot('S_total.CHG.dict.upstream.methylv.txt',ax2,10)
get_plot('S_total.CHH.dict.upstream.methylv.txt',ax3,10)
ax1.set_ylim(0,50)
ax2.set_ylim(0,50)
ax3.set_ylim(0,50)
plt.savefig('upstmet_exp_boxplot.png',dpi=300)


	



