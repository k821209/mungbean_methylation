#!/usr/bin/python3

file_in = 'total.fa.nonATGC.fa.out.gff'
dicT2info = {}
for line in open(file_in):
	cell = line.strip().split('\t') 
	strT = cell[-1]
	try:
		dicT2info[strT] += [line.strip()]
	except KeyError:
		dicT2info[strT] = [line.strip()]
for T in dicT2info:
	Outfile = open('TE_%s.gff'%T.replace('/','_'),'w')
	print('\n'.join(dicT2info[T]),file=Outfile)
	Outfile.close()

	


