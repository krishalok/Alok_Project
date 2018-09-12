from __future__ import print_function
import sys
import glob
import os

# path set to /projects/shilab/Hdesai/lung_cancer/RNASeq/


# os.chdir changes path to RNASeq
#os.chdir(path)


def create_header(files):
	#Will need to change IDs"
	header = "ID\t"
	if len(files) == 0:
		print ("No files found")
		sys.exit()
	
	#for filename in files:
        #sample_id = filename.split('-')[2]
        #header = header + sample_id + '\t'
    #header = header.rstrip()
    #return header

	for filename in files:
# 	sample_id = filename # changed to match file names 	
#		print(filename)
		sample_id = filename.split("TCGA")
		sample_id = sample_id[1]
		sample_id = sample_id[4:8]
#		print(sample_id)
# header = sample_id + "\t"
		header = header + sample_id + "\t"
#		print(header)
	header = header.rstrip()
	return header

def parse_genes(files):
    gene_data = {}
    gene_ids = []
    for data in files:
        genes = open(data, 'r').read()
        genes = genes.split('\n')
        for gene in genes:
            if not (gene.startswith('gene') or gene.startswith('?')) and gene != '':
                gene_id, value = gene.split('\t')
                gene_id = gene_id.split('|')[0]
                if gene_id in gene_data:
                    gene_data[gene_id] = gene_data[gene_id] + "\t" + value
                else:
                    gene_data[gene_id] = value
                    gene_ids.append(gene_id)

    return(gene_data, gene_ids)

#search argument needs to be quoted so it isn't expanded in the command line
def main():
    search = sys.argv[1]
    files = glob.glob(search)

    output = create_header(files)

    gene_data, gene_id = parse_genes(files)

    print(output)
    for gene in gene_id:
        print(gene + '\t' + gene_data[gene])

if __name__ == "__main__":
    main()
