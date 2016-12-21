#!/usr/bin/python3 

import sys,pickle
import numpy as np
file_in = sys.argv[1]
strChr = ''
def write_array2file(dic,Outfile_name): # arrayµçictµç
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()

def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)
i = 0
win = 1000
M = 0
U = 0
winML = {}
for line in open(file_in):
	cell = line.strip().split('\t')
	strChr = cell[0]
	strL	= cell[1]
	strM	= cell[3]
	strU	= cell[4]
	try:
		winML[(strChr,int(strL)/win)][1] += int(strU)
		winML[(strChr,int(strL)/win)][0] += int(strM)
	except KeyError:
		winML[(strChr,int(strL)/win)] = [int(strM),int(strU)]


MLs = np.array([M/(M+U+0.0001) for M,U in winML.values()])
		
write_array2file(MLs,file_in+'.wincount.npy')
		 		
