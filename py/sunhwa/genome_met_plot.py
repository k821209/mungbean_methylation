#!/usr/bin/python3 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
sns.set_palette("deep", desat=.6)
sns.set_context("poster",rc={"figure.figsize": (6,6)})
colors = sns.color_palette("deep", 3)


#/home/k821209/mungbean_methylation/sunhwa/S_1_4_5_1.fastq.gz_bismark_bt2_PE_report.txt
from matplotlib.gridspec import GridSpec
gs = GridSpec(3,2)
#labels = ['CpG', 'CHG', 'CHH', 'CN/CHN']
labels = ['CpG', 'CHG', 'CHH']

methylc_sunhwa = np.array([7804761, 9093042, 20106462]) # after pdf qv < 0.05
methylc_kyunggi = np.array([5924451,6612191,11916648]) # before binomial test
#unmethylc_sunhwa = np.array([92961756,146309079,1628774047,1055865]) # before binomial test
#unmethylc_sunhwa = np.array([4648741,7305309,70857661,1055865]) # before fdr pv < 0.001
unmethylc_sunhwa = np.array([5456665,8569864,91943365])
#unmethylc_kyunggi = np.array([51531259,80368953,902370054,1793849]) # before binomial test
methylc_sunhwa_chloro = np.array([76948,59426,57140])
unmethylc_sunhwa_chloro = np.array([24361518,20647334,127455973])

plt.figure(figsize=(12,12))

#explode = (0, 0, 0, 0.2) # only "explode" the 2nd slice (i.e. 'Hogs')
explode = (0, 0, 0)
plt.subplots_adjust(hspace = 0.4)

plt.subplot(gs[0,1])
plt.pie(methylc_sunhwa, explode=explode, labels=labels, colors=colors,
        autopct='%1.1f%%', shadow=True)
plt.title('S-Genome', fontsize=14, fontweight='bold')
plt.axis('equal')

plt.subplot(gs[1,1])
plt.pie(methylc_kyunggi, explode=explode, labels=labels, colors=colors,
        autopct='%1.1f%%', shadow=True)
plt.title('K-Genome', fontsize=14, fontweight='bold')
plt.axis('equal')


ax = plt.subplot(gs[0,0])
ind = np.arange(3)
width=0.35
bar1 = ax.bar(ind,methylc_sunhwa,width,color=colors[0])
bar2 = ax.bar(ind+width,unmethylc_sunhwa,width,color=colors[1])
ax.legend( (bar1[0], bar2[0]), ('MethylC', 'UnMethylC') )
plt.title('S-Genome', fontsize=14, fontweight='bold')
plt.ylabel('number of cytosine')
plt.ylim(0,2e8)
ax.set_xticks(ind+width)
ax.set_xticklabels(labels)
'''
ax = plt.subplot(gs[1,0])
ind = np.arange(4)
width=0.35
bar1 = ax.bar(ind,methylc_kyunggi,width,color='r')
bar2 = ax.bar(ind+width,unmethylc_kyunggi,width,color='b')
ax.legend( (bar1[0], bar2[0]), ('MethylC', 'UnMethylC') )
plt.title('K-Genome', fontsize=14, fontweight='bold')
plt.ylabel('number of cytosine')
plt.ylim(0,1e9)
ax.set_xticks(ind+width)
ax.set_xticklabels(labels)
'''
ax = plt.subplot(gs[1,0])
ind = np.arange(3)
width=0.35
bar1 = ax.bar(ind,methylc_sunhwa_chloro,width,color=colors[0])
bar2 = ax.bar(ind+width,unmethylc_sunhwa_chloro,width,color=colors[1])
ax.legend( (bar1[0], bar2[0]), ('MethylC', 'UnMethylC') )
plt.title('S-Chloroplast', fontsize=14, fontweight='bold')
plt.ylabel('number of cytosine')
plt.ylim(0,2e8)
ax.set_xticks(ind+width)
ax.set_xticklabels(labels)

plt.savefig("genome_met_plot.png",dpi=300)

plt.show()
