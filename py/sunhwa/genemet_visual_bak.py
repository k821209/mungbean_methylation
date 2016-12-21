#!/usr/bin/python3 
import sys,pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
sns.set_palette("deep", desat=.6)
sns.set_context("notebook",rc={"figure.figsize": (10, 2)})
colors = sns.color_palette("deep", 2)
def write_array2file(dic,Outfile_name):
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()

def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)

file_dict = sys.argv[1] # CG or CHG or CHH dict
win = 10000
dicLoc2Met = open_array(file_dict)

file_gff = 'Vradi_ver6.gff.sort.gff'

dicGN2loc = {}
dicID2GN = {}
for line in open(file_gff):
	if line.strip() == '':
		continue
	cell = line.strip().split('\t')
	strChr = cell[0]
	strL1 = cell[3]
	strL2 = cell[4]
	strT = cell[2]
	strInfo_list = cell[-1].rstrip(';').split(';')
	dicInfo = dict(zip([x.split('=')[0] for x in strInfo_list],[x.split('=')[1] for x in strInfo_list]))
	if strT == 'mRNA':
		strGN = dicInfo["Name"]
		strID = dicInfo["ID"]
		dicGN2loc[strGN] = (strChr,strL1,strL2,[])
		dicID2GN[strID] = strGN
	elif strT == 'exon':
		strID = dicInfo["Parent"]
		try:
			strGN = dicID2GN[strID]
		except KeyError:
			continue
		dicGN2loc[strGN][3].append([strL1,strL2])
target_gn = sys.argv[2]
Outfile_png = sys.argv[3] # name.png
strChr,strL1,strL2,cdslist = dicGN2loc[target_gn]
print(cdslist)
intL1 = int(strL1)
intL2 = int(strL2)
x = np.zeros(intL2-intL1+1)
indexlist1 = dicLoc2Met[(strChr,int(intL1/win))]
try:
	indexlist2 = dicLoc2Met[(strChr,int(intL1/win)-1)]
except KeyError:
	indexlist2 = []
try:
	indexlist3 = dicLoc2Met[(strChr,int(intL1/win)+1)]
except KeyError:
	indexlist3 = []
indexlist = indexlist1 + indexlist2 + indexlist3

for strLoc,strM,strU,strQ in indexlist:
	intL = int(strLoc)
	if intL1<=intL<=intL2:
		if (int(strM) + int(strU)) > 0:
			x[intL-intL1] = int(strM) / (int(strM) + int(strU))
#fig = plt.figure()
ax1 = plt.subplot()
ax1.plot(x,color=colors[0])
for l1,l2 in cdslist:
	intL1_cds = int(l1) - intL1
	intL2_cds = int(l2) - intL1
	ax1.add_patch(patches.Rectangle((intL1_cds, -0.03),intL2_cds-intL1_cds,0.03,color='r'))
plt.ylim(-0.03,1)
plt.title(target_gn)
plt.savefig(Outfile_png,dpi=300)
#plt.show()
