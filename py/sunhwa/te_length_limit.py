#!/usr/bin/python3 
import sys
file_gff = sys.argv[1] #'Vradi_ver6.gff.sort.gff'
Outfile = open(file_gff+'.genelimit1k.gff','w')
for line in open(file_gff):
	cell = line.strip().split('\t')
	if len(cell) < 3:
		continue
	strChr = cell[0]
	strL1 = cell[2]
	strL2 = cell[3]
	strT = cell[1]
	if strT == 'TE':
		if int(strL2) - int(strL1) > 1000:
			pass
		else:
			continue
		print(line.strip(),file=Outfile)	
