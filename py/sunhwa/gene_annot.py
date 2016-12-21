#!/usr/bin/python3 
import sys
file_annot= 'Vradi_ver6.annot'

dicGN2annt = {}
for line in open(file_annot):
	cell = line.strip().split('\t')
	strGN = cell[0]
	strAN = cell[-1]
	dicGN2annt[strGN] = strAN

for line in sys.stdin:
	cell = line.strip()
	strGN = cell
	try:
		print(strGN,dicGN2annt[strGN],sep='\t')
	except KeyError:
		print(strGN,'NA',sep='\t')
