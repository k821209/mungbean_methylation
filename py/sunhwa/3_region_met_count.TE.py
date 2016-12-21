#!/usr/bin/python3 
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import sys,pickle
import numpy as np
number_of_genes = int(sys.argv[3]) 
#number_of_genes = 1160
bin_number = 150 # genebody 5' 25 bins, 3' 25 bins, upstream 25bins, downstream 25bins
#bin_size = 100 # base pair  
genebody = np.zeros((number_of_genes,bin_number),dtype=float)
#genebody = np.zeros(bin_number)
#file_gff = 'Vradi_ver6.gff.sort.gff.test'
file_gff = sys.argv[2]#'Vradi_ver6.gff.sort.gff'
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

dicW2C = open_array(sys.argv[1]) # *.CG or *.CHG etc
i = 0
for line in open(file_gff):
	cell = line.strip().split('\t')
	if len(cell) < 3:
		continue
	strC = cell[0]
	strT = cell[1]
	#intL1 = int(cell[3]) - (bin_size*(bin_number/4))
	#intL2 = int(cell[4]) + (bin_size*(bin_number/4))
	intL1 = int(cell[2])
	intL2 = int(cell[3])
	if strT == 'TE':
#		print(i)
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
		pre_count = np.zeros(bin_number,dtype=float)
		bin_size = (intL2 - intL1) / 100
		if bin_size < 10:
			continue
		for j in range(25):
			M = 0
			U = 0
			for strLoc,strM,strU in indexed_list:
				if intL1-1000+(j*40) <= int(strLoc) < intL1-1000+ (j+1)*40:
					M += int(strM)
					U += int(strU)
				elif int(strLoc) > intL1:
					break
			T = M + U + 0.0001
			pre_count[j] = M/T
		for j in range(25):
			M = 0
			U = 0
			for strLoc,strM,strU in indexed_list:
				if intL2+(j*40) <= int(strLoc) < intL2+ (j+1)*40:
					M += int(strM)
					U += int(strU)
				elif int(strLoc) > intL2+1000:
					break
			T = M + U + 0.0001
			pre_count[j+125] = M/T
		for j in range(100):
			M = 0 
			U = 0
			for strLoc,strM,strU in indexed_list:
				if intL1+(j*bin_size) <= int(strLoc) < intL1+ (j+1)*bin_size:
					M += int(strM)
					U += int(strU)
				elif int(strLoc) > intL2:
					break
			#genebody[i,j] = count
			T = M + U + 0.0001
			pre_count[j+25] = M/T
#		print(T,pre_count,np.sum(pre_count))
		'''
		for j in range(int(bin_number/2)):
			M = 0 
			for strLoc,strM,strU in indexed_list:
				if intL2 - (j+1)*bin_size < int(strLoc) < intL2 - (j*bin_size):
					M += int(strM)
					U += int(strU)
				elif int(strLoc) > intL2:
					break
			#genebody[i,bin_number-1-j] = count
			pre_count[bin_number-1-j] = M/T
		'''
		#genebody[i] = pre_count/max(pre_count)
#		print(pre_count)
		genebody[i] = pre_count
		i += 1
#for each in genebody:
#	print(each)
#plt.boxplot(genebody, notch=0, sym='+', vert=1, whis=1.5)
write_array2file(genebody,sys.argv[1]+'.'+sys.argv[2]+'.100.npy')
'''
print(np.average(genebody,axis=1))
y = np.average(genebody.T,axis=1)
density = gaussian_kde(y)
density.covariance_factor = lambda : .25
density._compute_covariance()
x = np.arange(bin_number)
plt.plot(x,density(x))
plt.savefig(sys.argv[1] + '.test.png')
plt.show()	
'''
