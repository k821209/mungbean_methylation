#!/usr/bin/python3 
import sys
file_in = sys.argv[1]
for line in open(file_in):
	cell = line.strip().split()
	
