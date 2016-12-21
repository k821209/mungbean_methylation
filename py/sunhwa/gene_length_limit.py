#!/usr/bin/python3 

file_gff = 'Vradi_ver6.gff.sort.gff'
Outfile = open(file_gff+'.genelimit1k.gff','w')
for line in open(file_gff):
	cell = line.strip().split('\t')
	if len(cell) < 3:
		continue
	strChr = cell[0]
	strL1 = cell[3]
	strL2 = cell[4]
	strT = cell[2]
	if strT == 'gene':
		if int(strL2) - int(strL1) > 1000:
			pass
		else:
			continue
		print(line.strip(),file=Outfile)	
