#!/usr/bin/env python
# python 2.7.5 requires biopython

########### promer2circos ############


def rename_karyotype(ref, target, out_name):
    with open(ref, 'r') as f:
        contig_coords = []
        for row in f:
            print row
            data = row.rstrip().split(' ')
            if len(data) <3:
                data = row.rstrip().split('\t')
            contig_coords.append([data[3], int(data[4]), int(data[5])])
    data_list = read_circos_file(target)

    renamed_data = []
    for data in data_list:
        start = int(data[1])
        end = int(data[2])

        for one_contig in contig_coords:
            if start >= one_contig[1] and end <=one_contig[2]:
                data[0] = one_contig[0]
                renamed_data.append(data)
    with open(out_name, 'w') as new_circos_file:
        for row in renamed_data:
            new_circos_file.write('\t'.join(row)+'\n')


def read_circos_file(circos_file):
    data_list = []
    with open(circos_file) as f:
        for row in f:
            data = row.rstrip().split(' ')
            if len(data) <3:
                data = row.rstrip().split('\t')
            data_list.append(data)

    return data_list

if __name__ == '__main__':
    ###Argument handling.
    import argparse
    arg_parser = argparse.ArgumentParser(description='');
    #arg_parser.add_argument("coords_input", help="Directory to show-coords tab-delimited input file.");
    arg_parser.add_argument("-i", "--reference_karyotype", help="ref karyotype")
    arg_parser.add_argument("-t", "--target_karyotype", help="target karyotype")
    arg_parser.add_argument("-o", "--out_name", help="output name")

    args = arg_parser.parse_args()

    if not args.out_name:
        out_name = args.target_karyotype.split('.')[0] + '_renamed.' + args.target_karyotype.split('.')[1]

    rename_karyotype(args.reference_karyotype, args.target_karyotype, out_name)

