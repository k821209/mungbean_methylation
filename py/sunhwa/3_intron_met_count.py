#!/usr/bin/python3 
'''
Vr01    maker   mRNA    15780   22681   .       +       .       ID=329293;Name=Vradi01g00010.1;Parent=329292
Vr01    maker   gene    15780   22681   .       +       .       ID=329292;Name=Vradi01g00010;
Vr01    maker   five_prime_UTR  15780   16172   .       +       .       ID=329305;Name=Vradi01g00010.1;Parent=329293
Vr01    maker   exon    15780   16172   13921   +       .       ID=329294;Name=Vradi01g00010.1;Parent=329293
Vr01    maker   five_prime_UTR  16203   16296   .       +       .       ID=329306;Name=Vradi01g00010.1;Parent=329293
Vr01    maker   exon    16203   16296   13921   +       .       ID=329295;Name=Vradi01g00010.1;Parent=329293
Vr01    maker   five_prime_UTR  16568   16585   .       +       .       ID=329307;Name=Vradi01g00010.1;Parent=329293
Vr01    maker   exon    16568   16990   0.24    +       .       ID=329296;Name=Vradi01g00010.1;Parent=329293
Vr01    maker   start_codon     16586   16588   .       +       .       ID=329318;Name=Vradi01g00010.1;Parent=329293
Vr01    maker   CDS     16586   16990   .       +       0       ID=329308;Name=Vradi01g00010.1;Parent=329293
Vr01    maker   exon    17602   17662   0.24    +       .       ID=329297;Name=Vradi01g00010.1;Parent=329293
Vr01    maker   CDS     17602   17662   .       +       0       ID=329309;Name=Vradi01g00010.1;Parent=329293
Vr01    maker   exon    17772   17854   0.24    +       .       ID=329298;Name=Vradi01g00010.1;Parent=329293
'''
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import sys,pickle,kang
import numpy as np
win = 10000
def write_array2file(dic,Outfile_name): # arrayµçictµç
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()

def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)

file_fa = './gw/Vradi_ver6.fa'
dicHD2seq = kang.Fasta2dic(file_fa)

file_in = sys.argv[1]
file_gff = 'Vradi_ver6.gff.sort.gff.intronadded.gff'
file_fpkm = '/home/k821209/mungbean_methylation/RNAseq/diff_out/genes.fpkm_tracking' 
dicW2C = open_array(file_in) # *.CG or *.CHG etc

Outfile_kmer = open(file_in+'.cds.kmer.fa','w')
dkmer = {'CG':6,'CHG':7,'CHH':7}
i = 0
dicGN2count = {}
dicGN2fpkm = {}
for line in open(file_fpkm):
	cell = line.strip().split('\t')
	if cell[0] == 'tracking_id':
		continue
	strGN = cell[4].split(',')[-1]
	strFPKM_K = cell[9]
	strFPKM_S = cell[13]
	strFPKM_S_stat = cell[16]

	if strFPKM_S_stat == 'OK':
		pass
	else: continue
	try:
		if dicGN2fpkm[strGN]:
			print(1)
			exit()
	except KeyError:
		dicGN2fpkm[strGN] = strFPKM_S

dicID2mRNA = {} 

for line in open(file_gff):
	cell = line.strip().split('\t')
	if len(cell) < 3:
		continue
	strC = cell[0]
	strT = cell[2]
	strD = cell[6]
	intL1 = int(cell[3])
	intL2 = int(cell[4])
	if strT == 'mRNA':
		strInfo_list = cell[-1].rstrip(';').split(';')
		dicInfo = dict(zip([x.split('=')[0] for x in strInfo_list],[x.split('=')[1] for x in strInfo_list]))
		strID = dicInfo['ID']
		strGN = dicInfo['Name']
		dicID2mRNA[strID] = strGN
	if strT == 'intron':
		try:
			indexed_list1 = dicW2C[(strC,int(intL1/win)-1)]
		except KeyError:
			indexed_list1 = []
		try:
			indexed_list2 = dicW2C[(strC,int(intL1/win))]
		except KeyError:
			indexed_list2 = []
		try:
			indexed_list3 = dicW2C[(strC,int(intL1/win)+1)]
		except KeyError:
			indexed_list3 = []



		indexed_list = indexed_list1 + indexed_list2 + indexed_list3

		if indexed_list == []:
			i += 1
			continue
		strInfo_list = cell[-1].rstrip(';').split(';')
		dicInfo	= dict(zip([x.split('=')[0] for x in strInfo_list],[x.split('=')[1] for x in strInfo_list]))
		try:
			genename = dicID2mRNA[dicInfo['Parent']] 
		except KeyError:
			print(line)
			print(dicInfo['Parent'])
			exit()
		for strLoc,strM,strU,strDi in indexed_list:
			if intL1 <= int(strLoc) <= intL2:
				intLoc = int(strLoc)
				if sys.argv[2] == 'CG':
					if strDi == '+':
						context = dicHD2seq[strC][intLoc-1-2:intLoc+3]
					else:
						context = kang.rev_comp(dicHD2seq[strC][intLoc-1-3:intLoc+2])
				else:
					if strDi == '+':
						context = dicHD2seq[strC][intLoc-1-2:intLoc+4]
					else:
						context = kang.rev_comp(dicHD2seq[strC][intLoc-1-4:intLoc+2])
				if len(context) == dkmer[sys.argv[2]]:
					print('>'+strC+'_'+strLoc,file=Outfile_kmer)
					print(context,file=Outfile_kmer)
				else:
					pass
				try:
					dicGN2count[genename][1] += int(strU)
					dicGN2count[genename][0] += int(strM)
				except KeyError:
					dicGN2count[genename] = [int(strM),int(strU)]
Outfile = open(file_in +'.intronmethylv.txt','w')
print('genename','TotalM','TotalU','Methylation Level',sep='\t',file=Outfile)
for strGN in dicGN2count:
	intTotalM, intTotalU = dicGN2count[strGN]
	try:
		print(strGN, intTotalM, intTotalU, intTotalM/(intTotalM+intTotalU+0.00001),dicGN2fpkm[strGN.split('.')[0]],sep='\t',file=Outfile)
		
	except KeyError:
		pass
		#print(strGN, intTotalM, intTotalU, intTotalM/(intTotalM+intTotalU+0.00001),'NA',sep='\t',file=Outfile)

