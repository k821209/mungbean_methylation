#!/usr/bin/python3
#scaffold_7	106	+	0	13	CG	0
import pickle,sys

file_in = sys.argv[1]
dicW2C = {}
win = 10000
def write_array2file(dic,Outfile_name):
    Outfile = open(Outfile_name,'wb')
    pickle.dump(dic, Outfile)
    Outfile.close()

def open_array(filename):
    pkl_file = open(filename, 'rb')
    mydict2 = pickle.load(pkl_file)
    pkl_file.close()
    return(mydict2)

for line in open(file_in):
	cell = line.strip().split('\t')
	strC = cell[0]
	strL = cell[1]
	strD = cell[2]
	strM = cell[3]
	strU = cell[4]
	strQ = cell[-1]
	if 1:
		try:
			dicW2C[(strC,int(int(strL)/win))].append([strL,strM,strU,strD])
		except KeyError:
			dicW2C[(strC,int(int(strL)/win))] = [[strL,strM,strU,strD]]

write_array2file(dicW2C,file_in+'.q.dict')

