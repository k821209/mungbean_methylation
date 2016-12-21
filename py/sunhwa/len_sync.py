#!/usr/bin/python3 

import sys

seq_len = 7
file_in = sys.argv[1]
Outfile = open(file_in+'.lensync','w')
for line in open(file_in):
	seq = line.strip()
	if len(seq) == seq_len:
		print(line.strip(),file=Outfile)
