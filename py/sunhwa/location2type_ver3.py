#!/usr/bin/python3

import pickle,sys
'''
def write_dic2file(dic,Outfile_name):
        Outfile = open(Outfile_name,'wb')
        pickle.dump(dic, Outfile)
        Outfile.close()
def open_dicfile(filename):
        pkl_file = open(filename, 'rb')
        mydict2 = pickle.load(pkl_file)
	strChr	= cell[0]
        pkl_file.close()

        return(mydict2)
'''
#Vr01	306	306	CG	0.909090909091	10	1
file_in		= sys.argv[1] #'all.out.soja.mf.polymorphic' 
file_gff 	= 'Vradi_ver6.gff.sort.gff'
Outfile		= open(file_in+'.loc2tp','w')
Outfile_mm	= open('multiple_match.txt','w')

intW	= 100000
dicCL2line	= {} # Key = (strChr, intLoc), Value = line
i = 0
dicPACID2TGN = {}
for line in open(file_gff):
	if line[0] == '#' or line.strip() == '':
		continue
	cell = line.strip().split('\t')
	strChr	= cell[0]
	strLoc1	= cell[3]
	intLoc	= int(int(strLoc1)/intW)
	strTP = cell[2] 
	if strTP == 'mRNA':
		strInfo_list = cell[-1].rstrip(';').split(';')
		dicInfo = dict(zip([x.split('=')[0] for x in strInfo_list],[x.split('=')[1] for x in strInfo_list]))
		pacid = dicInfo['ID']
		name = dicInfo['Name']
		dicPACID2TGN[pacid] = name
	if i % 100000 == 0:
		print('processing : %d'%i,end='\r')
	i += 1
	try:
		dicCL2line[(strChr,intLoc)].append(line.strip())
	except KeyError:
		dicCL2line[(strChr,intLoc)] = [line.strip()]      # (strChr, intLoc) and (strChr, intLoc-1) should be merged for genes ranging like 199999(1) - 200011(2) harboring SNP 200000(2)
#write_dic2file(dicCL2line,file_gff_bdic)

for line in open(file_in):
	if line[0] == '#':
		continue
	cell = line.strip().split('\t')
	strChr = cell[0]
	strLoc1 = cell[1]
	#print(strChr,strLoc1)
	intLoc	= int(int(strLoc1)/intW)
	if intLoc-1 < 0:
		try:
			search_list 	= dicCL2line[(strChr,intLoc)]
		except KeyError:
			print(line.strip(),"Intergenic",'NA',sep='\t',file=Outfile)
			continue
	else:	
		try:
			search_list	= dicCL2line[(strChr,intLoc-1)] + dicCL2line[(strChr,intLoc)]
		except KeyError:
			try:
				search_list	= dicCL2line[(strChr,intLoc-1)]
			except KeyError:
				try:
					search_list     = dicCL2line[(strChr,intLoc)]
				except KeyError:
					print(line.strip(),"Intergenic",'NA',sep='\t',file=Outfile)	
					continue	
				
	bmRNA	= 0
	bCDS	= 0
	b1kup	= 0
	b5UTR	= 0
	b3UTR	= 0
	dicGN2tp = {}
	for each in search_list:
		each_cell 	= each.split('\t')
		strTchr		= each_cell[0]
		strTstrand	= each_cell[6]
		strTtyp		= each_cell[2]
		strTloc1	= each_cell[3]
		strTloc2	= each_cell[4]
		strInfo_list = each_cell[-1].rstrip(';').split(';')
		dicInfo = dict(zip([x.split('=')[0] for x in strInfo_list],[x.split('=')[1] for x in strInfo_list]))
		if strTtyp == 'gene' or strTtyp == 'contig':
			continue
		elif strTtyp == 'mRNA':
			strTGN		= dicInfo['Name']
		else :
			try:
				strTGN		= dicPACID2TGN[dicInfo['Parent']]
			except KeyError:
				print(each)
				exit()
	#	print(strTGN)	
		if strTGN[-2:] != '.1':
			continue
		if strTtyp == 'mRNA':
			strTGN          = dicInfo['Name']
			if strTstrand == '+' and int(strTloc1)-1000 <= int(strLoc1) <= int(strTloc1):
				b1kup 	= 1
				try:
					dicGN2tp[strTGN].append('1kup')		
				except KeyError:
					dicGN2tp[strTGN] = ['1kup']
			if strTstrand == '-' and int(strTloc2) <= int(strLoc1) <= int(strTloc2) + 1000:
				b1kup   = 1
				try:
					dicGN2tp[strTGN].append('1kup')		
				except KeyError:
					dicGN2tp[strTGN] = ['1kup']
				
		if int(strTloc1) <= int(strLoc1) <= int(strTloc2):
			if strTtyp == 'mRNA':
				bmRNA 	= 1
				try:
					dicGN2tp[strTGN].append('mRNA')		
				except KeyError:
					dicGN2tp[strTGN] = ['mRNA']
			elif strTtyp == 'CDS':
				bCDS 	= 1
				try:
					dicGN2tp[strTGN].append('CDS')		
				except KeyError:
					dicGN2tp[strTGN] = ['CDS']
			elif strTtyp == 'five_prime_UTR':
				b5UTR	= 1
				try:
					dicGN2tp[strTGN].append('5UTR')		
				except KeyError:
					dicGN2tp[strTGN] = ['5UTR']
			elif strTtyp == 'three_prime_UTR':
				b3UTR	= 1
				try:
					dicGN2tp[strTGN].append('3UTR')		
				except KeyError:
					dicGN2tp[strTGN] = ['3UTR']
		if int(strLoc1)+10000 < int(strTloc1):
			break
	dicGN2tp_list = list(dicGN2tp.keys())
	if len(dicGN2tp_list) > 1:
		print(line.strip(),dicGN2tp,file=Outfile_mm)
		continue
	if len(dicGN2tp_list) == 0:
		strMatched_GN = 'NA'
	else: strMatched_GN = dicGN2tp_list[0]
	if bmRNA == 1 and bCDS == 1:
		#print(line.strip(),"CDS",strMatched_GN,sep='\t')
		print(line.strip(),"CDS",strMatched_GN,sep='\t',file=Outfile)
	elif bmRNA == 0 and b1kup == 1:
		#print(line.strip(),"1_kb_upstream",strMatched_GN,sep='\t')
		print(line.strip(),"1_kb_upstream",strMatched_GN,sep='\t',file=Outfile)
	elif bmRNA == 1 and b5UTR == 1:
		#print(line.strip(),"five_prime_UTR",strMatched_GN,sep='\t')
		print(line.strip(),"five_prime_UTR",strMatched_GN,sep='\t',file=Outfile)
	elif bmRNA == 1 and b3UTR == 1:
		#print(line.strip(),"three_prime_UTR",strMatched_GN,sep='\t')
		print(line.strip(),"three_prime_UTR",strMatched_GN,sep='\t',file=Outfile)
	elif bmRNA == 1:
		#print(line.strip(),"Intron",strMatched_GN,sep='\t')
		print(line.strip(),"Intron",strMatched_GN,sep='\t',file=Outfile)
	else: 
		#print(line.strip(),"Intergenic",strMatched_GN,sep='\t')
		print(line.strip(),"Intergenic",dicGN2tp,sep='\t',file=Outfile)


	
	
	

