#!/usr/bin/python3 
import sys
file_annot = '../Vradi_ver6.annot'
dicGN2GO = {}
for line in open(file_annot):
	if line.strip() == '':
		continue	
	cell = line.strip().split('\t')
	if len(cell) < 5:
		continue
	try:
		strGN = cell[0]
		strGO = cell[4]
	except :
		print(line)
		exit()
	
	
	if strGO == '':
		continue
	#print(strGO)
	dicGN2GO[strGN] = strGO

for line in sys.stdin:
	strGN = line.strip()
	try:
		print('\n'.join(dicGN2GO[strGN].split(',')))
	except KeyError:
		pass
