# This script splits the unified gtex rpkm data matrix by tissue type

import pandas as pd
import sys

fname_rna = sys.argv[1] # unified RNAseq file (gct file)
fname_sampleann = sys.argv[2] # sample annotation file
outfolder = sys.argv[3] # save files in this folder

rna_exp = pd.read_csv(fname_rna,header=2,sep='\t')
sample_ann = pd.read_csv(fname_sampleann,sep='\t',header=0)

tissue_list = sample_ann.ix[:,1]
tissue_list = tissue_list.drop_duplicates()
tissue_list = tissue_list.tolist()
allsamples = rna_exp.columns.values

# split for each unique tissue
for i in tissue_list:
	current_tissue = i
	print current_tissue
	current_samples = sample_ann[sample_ann.SMTSD==current_tissue]
	current_samples = current_samples.ix[:,0].tolist()
	current_samples = list(set(current_samples).intersection(allsamples))
	current_rna_exp = rna_exp[current_samples]
	current_rna_exp['Name']=rna_exp.ix[:,0]
	current_rna_exp['Description']=rna_exp.ix[:,1]
	r,n = current_rna_exp.shape
	current_rna_exp = current_rna_exp.ix[:,[n-2,n-1]+range(n-2)]
	fname = current_tissue.replace(" ","")+'.txt'
	fname = outfolder+fname
	current_rna_exp.to_csv(fname,sep='\t',index=False)
	del current_rna_exp, current_samples
