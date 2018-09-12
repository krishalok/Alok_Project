from __future__ import print_function
import sys
import glob

#TODO: Will have to check TCGA barcoding standard for 01A vs 01B vs 11A etc

def create_header(files):
    header = 'ID\t'

    if len(files) == 0:
        print('No files found in the glob')
        sys.exit()

    for filename in files:
        sample_id = filename.split('-')[2]
        header = header + sample_id + '\t'
    header = header.rstrip()
    return header

def parse_mirna(files):
    mirna_data = {}
    mirna_ids = []

    for counter, filename in enumerate(files):
        mirna_file = open(filename, 'r').read()
        mirnas = mirna_file.split('\n')
        for mirna in mirnas:
            if (not mirna.startswith('miR')) and mirna != '':
                mirna_id, _, value, _ = mirna.split('\t')
                if mirna_id in mirna_data:
                    data_len = len(mirna_data[mirna_id].split('\t'))
                    if (data_len != (counter)):
                        mirna_data[mirna_id] = mirna_data[mirna_id] + '\tNA'*(counter-data_len) + '\t' + value
                    else:
                        mirna_data[mirna_id] = mirna_data[mirna_id] + '\t' + value
                else:
                    if counter==1:
                        mirna_data[mirna_id] = value
                    else:
                        mirna_data[mirna_id] = 'NA\t'*(counter) + value
                    mirna_ids.append(mirna_id)

    for mirna in mirna_ids:
        mirna_line = mirna_data[mirna]
        expected_length = len(files)
        actual_length = len(mirna_line.split('\t'))
        if actual_length != expected_length:
            diff = expected_length - actual_length
            mirna_data[mirna] = mirna_line + '\tNA'*diff


    return(mirna_data, mirna_ids)

def main():
    search = sys.argv[1]
    files = glob.glob(search)
    files = [file for file in files if '-01A' in file]

    output = create_header(files)

    mirna_data, mirna_ids = parse_mirna(files)

    print(output)
    for mirna in mirna_ids:
        print(mirna + '\t' + mirna_data[mirna])

if __name__ == '__main__':
    main()
