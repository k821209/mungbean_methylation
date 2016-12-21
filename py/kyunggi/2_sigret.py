#!/usr/bin/python3
#Vr01    5       -       0       0       CHH     1.0     1 1
import sys

file_in = sys.argv[1] #*.qv
Outfile = open(file_in+'.qvsig.txt','w')
for line in open(file_in):
	cell = line.strip().split()
	strS,strL,strD,strM,strU,strT,strP,strP1,strQ = cell
	if float(cell[-1]) < 0.05:
		print('\t'.join(cell),file=Outfile)
	else : print(strS,strL,strD,0,int(strM) + int(strU),strT,strP,strP1,strQ,sep='\t',file=Outfile)
