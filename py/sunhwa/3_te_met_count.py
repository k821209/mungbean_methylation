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
file_gff = sys.argv[2]
dicW2C = open_array(file_in) # *.CG or *.CHG etc
i = 0
dicGN2count = {}
for line in open(file_gff):
	cell = line.strip().split('\t')
	if len(cell) < 3:
		continue
	strC = cell[0]
	strT = cell[1]
	strD = '+'
	intL1 = int(cell[2])
	intL2 = int(cell[3])
	if strT == 'TE':
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
		genename = cell[-1] + '.' + strC + '_' + str(intL1)
		for strLoc,strM,strU in indexed_list:
			if intL1 <= int(strLoc) <= intL2:
				try:
					dicGN2count[genename][1] += int(strU)
					dicGN2count[genename][0] += int(strM)
				except KeyError:
					dicGN2count[genename] = [int(strM),int(strU)]
Outfile = open(file_in + '.'+ file_gff +'.genemethylv.txt','w')
print('genename','TotalM','TotalU','Methylation Level',sep='\t',file=Outfile)
for strGN in dicGN2count:
	intTotalM, intTotalU = dicGN2count[strGN]
	try:
		print(strGN, intTotalM, intTotalU, intTotalM/(intTotalM+intTotalU+0.00001),sep='\t',file=Outfile)
		
	except KeyError:
		pass
		#print(strGN, intTotalM, intTotalU, intTotalM/(intTotalM+intTotalU+0.00001),'NA',sep='\t',file=Outfile)

