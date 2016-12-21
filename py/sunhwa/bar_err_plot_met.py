#!/usr/bin/python3 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
sns.set_palette("deep")
sns.set_context("poster",rc={"figure.figsize": (8,6)})
boxColors = sns.color_palette("deep", 3)

file_gm_CG = 'S_total.CG.q.dict.exonmethylv.txt'
file_gm_CHG = 'S_total.CHG.q.dict.exonmethylv.txt'
file_gm_CHH = 'S_total.CHH.q.dict.exonmethylv.txt'
data_CG = np.loadtxt(file_gm_CG,usecols = (3,),skiprows=1)
data_CHG = np.loadtxt(file_gm_CHG,usecols = (3,),skiprows=1)
data_CHH = np.loadtxt(file_gm_CHH,usecols = (3,),skiprows=1)

data_CG = data_CG[(data_CG>0)]
data_CHG = data_CHG[(data_CHG>0)]
data_CHH = data_CHH[(data_CHH>0)]

sns.violinplot([data_CG,data_CHG,data_CHH],names=["CG","CHG","CHH"],widths=0.8)
#bp = plt.boxplot([data_CG,data_CHG,data_CHH], notch=0,showfliers=False, sym='+', vert=1, whis=1.5)
#plt.setp(bp['boxes'], color='black')
#plt.setp(bp['whiskers'], color='black')
#plt.setp(bp['fliers'], color='red', marker='+')


plt.ylabel('Weighted methylation level')
plt.savefig('exonmet.violinplot.png',dpi=300)


