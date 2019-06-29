#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------

def gbk2fna(seq_records, outname=False):
    from Bio.SeqRecord import SeqRecord

    if type(seq_records) == str:
        record_list = [i for i in SeqIO.parse(open(seq_records), 'genbank')]
    elif type(seq_records) == list and type(seq_records[0])==str:
        record_list = []
        for i in seq_records:
            tmp_list = [i for i in SeqIO.parse(open(i), 'genbank')]
            record_list+=tmp_list

    elif isinstance(seq_records, SeqRecord):
        record_list = [seq_records]
    elif type(seq_records) == list and isinstance(seq_records[0], SeqRecord):
        record_list = seq_records
    else:
        print ('wrong inpur reference')

    length_records = [len(i.seq) for i in record_list]
    longest_record = length_records.index(max(length_records))

    if not outname:
        outname = record_list[longest_record].id.split('.')[0] + ".fna"

    max_len = 0
    if len(record_list) > 1:
        rec_list = []
        output_handle = open(outname, "w")
        for record in record_list:
            output_handle.write(">%s %s\n%s\n" % (record.name, record.description,record.seq))
    else:
        output_handle = open(outname, "w")
        try:
            gi = record_list[0].annotations["gi"]
            source = record_list[0].annotations["source"]
            output_handle.write(">%s gi|%s|%s\n%s\n" % (record_list[0].name, gi, source, record_list[0].seq))
        except:
            output_handle.write(">%s\n%s\n" % (record_list[0].name, record_list[0].seq))

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)


    args = parser.parse_args()
    input_handle = open(args.input_gbk, "rU")
    seq_records = list(SeqIO.parse(input_handle, "genbank"))

    gbk2fna(seq_records, args.outname)
