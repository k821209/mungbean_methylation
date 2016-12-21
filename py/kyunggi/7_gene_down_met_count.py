#!/usr/bin/python3 
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import sys,pickle
import numpy as np
win = 10000
def write_array2file(dic,Outfile_name): # arrayµçictµç
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()

def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)

file_in = sys.argv[1]
file_gff = 'Vradi_ver6.gff.sort.gff'
file_fpkm = '/home/k821209/mungbean_methylation/RNAseq/diff_out/genes.fpkm_tracking' 
dicW2C = open_array(file_in) # *.CG or *.CHG etc
i = 0
dicGN2count = {}
dicGN2fpkm = {}
for line in open(file_fpkm):
	cell = line.strip().split('\t')
	if cell[0] == 'tracking_id':
		continue
	strGN = cell[4].split(',')[-1]
	strFPKM_K = cell[10]
	strFPKM_S = cell[14]
	strFPKM_S_stat = cell[16]

	if strFPKM_S_stat == 'OK':
		pass
	else: continue
	try:
		if dicGN2fpkm[strGN]:
			print(1)
			exit()
	except KeyError:
		dicGN2fpkm[strGN] = strFPKM_K



for line in open(file_gff):
	cell = line.strip().split('\t')
	if len(cell) < 3:
		continue
	strC = cell[0]
	strT = cell[2]
	strD = cell[6]
	intL1 = int(cell[3])
	intL2 = int(cell[4])
	if strT == 'gene':
		try:
			indexed_list1 = dicW2C[(strC,int(intL1/win)-1)]
		except KeyError:
			indexed_list1 = []
		try:
			indexed_list2 = dicW2C[(strC,int(intL1/win))]
		except KeyError:
			indexed_list2 = []
		try:
			indexed_list3 = dicW2C[(strC,int(intL1/win)+1)]
		except KeyError:
			indexed_list3 = []



		indexed_list = indexed_list1 + indexed_list2 + indexed_list3

		if indexed_list == []:
			i += 1
			continue
		strInfo_list = cell[-1].split(';')[:-1]
		dicInfo	= dict(zip([x.split('=')[0] for x in strInfo_list],[x.split('=')[1] for x in strInfo_list]))
		genename = dicInfo['Name'] 
		for strLoc,strM,strU,strQ in indexed_list:
			if strD == '-':
				if intL1-1000 <= int(strLoc) <= intL1:
					try:
						dicGN2count[genename][1] += int(strU)
						dicGN2count[genename][0] += int(strM)
					except KeyError:
						dicGN2count[genename] = [int(strM),int(strU)]
			else: 
				if intL2 <= int(strLoc) <= intL2+1000:
					try:
						dicGN2count[genename][1] += int(strU)
						dicGN2count[genename][0] += int(strM)
					except KeyError:
						dicGN2count[genename] = [int(strM),int(strU)]
					
Outfile = open(file_in +'.downstream.methylv.txt','w')
print('genename','TotalM','TotalU','Methylation Level',sep='\t',file=Outfile)
for strGN in dicGN2count:
	intTotalM, intTotalU = dicGN2count[strGN]
	try:
		print(strGN, intTotalM, intTotalU, intTotalM/(intTotalM+intTotalU+0.00001),dicGN2fpkm[strGN],sep='\t',file=Outfile)
		
	except KeyError:
		pass
		#print(strGN, intTotalM, intTotalU, intTotalM/(intTotalM+intTotalU+0.00001),'NA',sep='\t',file=Outfile)

