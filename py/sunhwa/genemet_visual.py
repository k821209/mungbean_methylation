#!/usr/bin/python3 
import sys,pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
sns.set_palette("deep", desat=.6)
sns.set_context("poster",rc={"figure.figsize": (10, 4)})
colors = sns.color_palette("deep", 3)
def write_array2file(dic,Outfile_name):
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()

def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)

file_dict = sys.argv[1] # CG or CHG or CHH dict header
win = 10000
dicLoc2MetCG = open_array(file_dict+'.CG.q.dict')
dicLoc2MetCHG = open_array(file_dict+'.CHG.q.dict')
dicLoc2MetCHH = open_array(file_dict+'.CHH.q.dict')

file_gff = 'Vradi_ver6.gff.sort.gff'

dicGN2loc = {}
dicID2GN = {}
for line in open(file_gff):
	if line.strip() == '':
		continue
	cell = line.strip().split('\t')
	if len(cell) < 3:
		continue
	strChr = cell[0]
	strL1 = cell[3]
	strL2 = cell[4]
	strD = cell[6]
	strT = cell[2]
	strInfo_list = cell[-1].rstrip(';').split(';')
	dicInfo = dict(zip([x.split('=')[0] for x in strInfo_list],[x.split('=')[1] for x in strInfo_list]))
	if strT == 'mRNA':
		strGN = dicInfo["Name"]
		strID = dicInfo["ID"]
		dicGN2loc[strGN] = (strChr,strL1,strL2,strD,[])
		dicID2GN[strID] = strGN
	elif strT == 'exon':
		strID = dicInfo["Parent"]
		try:
			strGN = dicID2GN[strID]
		except KeyError:
			continue
		dicGN2loc[strGN][4].append([strL1,strL2])

def get_x(dicLoc2Met,intL1,intL2):
	x = np.zeros(intL2-intL1+1+2000)
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

	for strLoc,strM,strU,strD_met in indexlist:
		intL = int(strLoc)
		if intL1-1000<=intL<=intL2+1000:
			if (int(strM) + int(strU)) > 0:
				if strD_met == '+':
					x[intL-intL1+1000] = int(strM) / (int(strM) + int(strU))
				elif strD_met == '-':
					x[intL-intL1+1000] = -1 * (int(strM) / (int(strM) + int(strU)))
				else:
					print(1)
					exit()
	if strD == '+':
		pass
	else:
		x = x[::-1] * (-1)
	return(x)

target_gn_list = sys.argv[2].split(',')
header = sys.argv[3] 
for target_gn in target_gn_list:
	strChr,strL1,strL2,strD,cdslist = dicGN2loc[target_gn+'.1']
	print(cdslist)
	intL1 = int(strL1)
	intL2 = int(strL2)
	x_CG = get_x(dicLoc2MetCG,intL1,intL2)
	x_CHG = get_x(dicLoc2MetCHG,intL1,intL2)
	x_CHH = get_x(dicLoc2MetCHH,intL1,intL2)


	fig, (ax3,ax2,ax1) = plt.subplots(3,1,sharex=True)

	#ax1 = fig.add_subplot(gs[2,0])
	#ax2 = fig.add_subplot(gs[1,0],sharex=ax1)
	#ax3 = fig.add_subplot(gs[0,0],sharex=ax1)
	ind = np.arange(len(x_CG))
	rects1 = ax3.bar(ind,x_CG,1,
                color=colors[0],edgecolor='none') 
	rects2 = ax2.bar(ind,x_CHG,1,
                color=colors[1],edgecolor='none') 
	rects3 = ax1.bar(ind,x_CHH,1,
                color=colors[2],edgecolor='none') 


	#ax1.plot(x,color=colors[0])
	genelength = len(x_CG)
	for l1,l2 in cdslist:
		intL1_cds = int(l1) - intL1
		intL2_cds = int(l2) - intL1
		if strD == '+':
			ax1.add_patch(patches.Rectangle((intL1_cds+1000, -0.04),intL2_cds-intL1_cds,0.03,color='r'))
			ax2.add_patch(patches.Rectangle((intL1_cds+1000, -0.04),intL2_cds-intL1_cds,0.03,color='r'))
			ax3.add_patch(patches.Rectangle((intL1_cds+1000, -0.04),intL2_cds-intL1_cds,0.03,color='r'))
		elif strD == '-':
			ax1.add_patch(patches.Rectangle((genelength-(intL2_cds+1000), -0.04),intL2_cds-intL1_cds,0.03,color='r'))
			ax2.add_patch(patches.Rectangle((genelength-(intL2_cds+1000), -0.04),intL2_cds-intL1_cds,0.03,color='r'))
			ax3.add_patch(patches.Rectangle((genelength-(intL2_cds+1000), -0.04),intL2_cds-intL1_cds,0.03,color='r'))
		else:
			print(1)
			exit()
	ax1.set_ylim(-1,1)
	ax2.set_ylim(-1,1)
	ax3.set_ylim(-1,1)
	ax3.set_title(target_gn)
	#ax3.legend((rects1[0], rects2[0], rects3[0]), ('CG', 'CHG','CHH') )
	plt.savefig(header + '.' + target_gn+'.png',dpi=300)
#plt.show()
