#!/usr/bin/python3 

import kang,sys

dicHD2seq = kang.Fasta2dic(sys.argv[1])

print(dicHD2seq[sys.argv[2]][int(sys.argv[3])-1:int(sys.argv[4])])
