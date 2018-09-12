import os, sys, re

new_dict = {}

input_file = open('/projects/shilab/Data/TCGA/COCA/Clinical/aliquot.txt')
for line in input_file:
#	print "these are the lines:", line
	new_name, old = line.split()
	new_dict[old] = new_name

#print new_dict
#os.chdir("/projects/shilab/Hdesai/lung_cancer/RNASeq/")

path = "/projects/shilab/Data/TCGA/COCA/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/test/"
os.chdir(path)
for filename in os.listdir(path):
	root, extension = os.path.splitext(filename)
	extra = root[36:]
	file_end = root[:36]
#	print "this is the extension", extension
	if file_end in new_dict:
		print "this works!"
# old name contains extra and keep --- keep needs to be kept for renaming	
		old = path + filename
#		print "this is the old name:", old
		new = path+new_dict[file_end]+ extra+extension
#		print "this will be the new name", new
	#	print old, new
		os.rename(old, new)
