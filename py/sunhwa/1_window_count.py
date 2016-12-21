#!/usr/bin/python3 
#scaffold_7	310	+	2	9	CHH	6.7826031573e-05	6.782603e-05	0.0003467132
import sys

file_in = sys.argv[1]

window = 50000 # 100Kb window
dicW2pC = {}
dicW2mC = {}
for line in open(file_in):
	strChr, strL, strD, strM, strU, strT, strP, strP1, strQ = line.strip().split('\t')
	if float(strQ) < 0.05:
		strMet = 1
	else:
		strMet = 0

	if strD == '+':
		try:
			dicW2pC[int(int(strL)/window)] += int(strMet)
		except KeyError:
			dicW2pC[int(int(strL)/window)] = int(strMet)
	else: 
		try:
			dicW2mC[int(int(strL)/window)] += int(strMet) * -1
		except KeyError:
			dicW2mC[int(int(strL)/window)] = int(strMet) * -1
dicW2pC_list = list(dicW2pC.keys())
dicW2pC_list.sort()
dicW2mC_list = list(dicW2mC.keys())
dicW2mC_list.sort()
Outfile1 = open(file_in+'.phist','w')
Outfile2 = open(file_in+'.mhist','w')
for intW in dicW2pC_list:
	print(intW,dicW2pC[intW],sep='\t',file=Outfile1)
for intW in dicW2mC_list:
	print(intW,dicW2mC[intW],sep='\t',file=Outfile2)

