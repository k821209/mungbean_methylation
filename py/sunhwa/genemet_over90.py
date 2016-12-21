#!/usr/bin/python3
#S_total.CG.q.dict.genemethylv.txt
import sys

file_in = 'S_total.CG.q.dict.genemethylv.txt'
for line in open(file_in):
	cell = line.strip().split('\t')
	if line[0:8] == 'genename':
		continue
	strGN = cell[0]
	strM  = cell[1]
	strU  = cell[2]
	strML = cell[3]
	if float(strML) > 0.90:
		print(strGN + '.1')	
