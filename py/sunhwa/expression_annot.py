#!/usr/bin/python3

import sys

dicGN2ex = {}
dicGN2ml = {}
file_ml = sys.argv[1]#'S_total.CHG.dict.genemethylv.txt'
for line in open(file_ml):
	cell = line.strip().split('\t')
	if cell[0] == 'genename':
		continue
	strGN = cell[0]
	strML = cell[3]
	dicGN2ml[strGN] = strML
file_ex = '/home/k821209/mungbean_methylation/RNAseq/diff_out/genes.fpkm_tracking'
for line in open(file_ex):
	cell = line.strip().split('\t')
	if cell[0] == 'tracking_id':
		continue
	strGN = cell[4].split(',')[-1]
	strEXP = cell[13]
	if cell[-1] == 'OK':
		dicGN2ex[strGN] = strEXP
for line in sys.stdin:
	strGN = line.strip()
	try:
		strML = dicGN2ml[strGN]
	except KeyError:
		strML = 'NA'
	try:
		print(strGN,strML,dicGN2ex[strGN],sep='\t')
	except KeyError:
		print(strGN,strML,'NA',sep='\t')
