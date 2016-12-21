'''
 CG
 154728 1_kb_upstream
  86523 CDS
   4120 five_prime_UTR
7206642 Intergenic
 341803 Intron
   3406 three_prime_UTR
 CHG
 149127 1_kb_upstream
  83158 CDS
   4167 five_prime_UTR
8449689 Intergenic
 395491 Intron
   3768 three_prime_UTR
 CHH
 462398 1_kb_upstream
 116722 CDS
   8205 five_prime_UTR
18677479 Intergenic
 819679 Intron
   7379 three_prime_UTR
'''
#!/usr/bin/python3 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
sns.set_palette("deep", desat=.6)
sns.set_context("poster",rc={"figure.figsize": (10,10)})

colors = sns.color_palette("deep", 3)

x = np.arange(6)
CG = np.array([154728,86523,4120,7206642,341803,3406])
CHG = np.array([149127,83158,4167,8449689,395491,3768])
CHH = np.array([462398,116722,8205,18677479,819679,7379])
CG = CG/np.sum(CG)
CHG = CHG/np.sum(CHG)
CHH = CHH/np.sum(CHH)
width = 0.3

ax = plt.subplot(111)
ax.bar(x,CG,color=colors[0],width=0.3,label='CG')
ax.bar(x+width,CHG,color=colors[1],width=0.3,label='CHG')
ax.bar(x+width+width,CHH,color=colors[2],width=0.3,label='CHH')
ax.set_xticks(x+width)
ax.set_xticklabels(['Upstream','CDS','5`UTR','Intergenic','Intron','3`UTR'],rotation=40)

#plt.show()
plt.legend()
plt.savefig('classification.met.png',dpi=300)

