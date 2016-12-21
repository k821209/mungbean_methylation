#!/usr/bin/python3 
#scaffold_7      106     +       0       13      CG      1.0     1       1
import kang,sys

file_fa = './gw/Vradi_ver6.fa'
file_in = sys.argv[1]
dicHD2seq = kang.Fasta2dic(file_fa)
Outfile = open(file_in+'.context','w')
Outfile_fa = open(file_in+'.context.fa','w')
seq_len = 6
for line in open(file_in):
	cell = line.strip().split('\t')
	strChr = cell[0]

	intLoc = int(cell[1])
	strD = cell[2]
	if strD == '+':
		context = dicHD2seq[strChr][intLoc-1-2:intLoc+3]
	else:
		context = kang.rev_comp(dicHD2seq[strChr][intLoc-1-3:intLoc+2])
	if len(context) == seq_len:
		print('>'+strChr+'_'+cell[1],file=Outfile_fa)
		print(context,file=Outfile_fa)
	print(line.strip(),context,sep='\t',file=Outfile)
