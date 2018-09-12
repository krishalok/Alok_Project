from __future__ import print_function
import sys
import glob

def create_header(files):
    #Will need to change IDs
    header = 'ID\t'

    if len(files) == 0:
        print('No files found in the glob')
        sys.exit()

    for filename in files:
        sample_id = filename.split('-')[2]
        header = header + sample_id + '\t'
    header = header.rstrip()
    return header

def parse_genes(files):
    gene_data = {}
    gene_ids = []
    for counter, data in enumerate(files):
        genes = open(data, 'r').read()
        genes = genes.split('\n')
        for gene in genes:
            if not (gene.startswith('gene') or gene.startswith('?')) and gene != '':
                gene_id, value = gene.split('\t')
                orig_id = gene_id
                gene_id = gene_id.split('|')[0]
                if orig_id in gene_data:
                    gene_data[orig_id] = gene_data[orig_id] + '\t' + value
                elif gene_id in gene_data and counter!=0:
                    gene_data[gene_id] = gene_data[gene_id] + '\t' + value
                elif gene_id in gene_data and counter==0:
                    gene_data[orig_id] = value
                    gene_ids.append(orig_id)
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

if __name__ == '__main__':
    main()
