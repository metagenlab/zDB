#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert gbk file to ffn
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------



def gbk2ffn(seq_records, outname, locus_tag=False, genome_accession=False):
    from Bio.SeqRecord import SeqRecord

    if type(seq_records) == list:
        if type(seq_records[0]) == str:
           seq_records = list(SeqIO.parse(input_handle, "genbank"))
        elif isinstance(seq_records[0], SeqRecord):
            pass
    elif isinstance(seq_records, SeqRecord):
        seq_records = [seq_records]
    else:
        raise('unknown input format of record')

    max_len = 0
    if len(seq_records) > 1:
        rec_list = []
        for record in seq_records:
            if len(record.seq)>max_len:
                print('plus long!!!!!!!!!!!!!!')
                max_len = len(record.seq)
                print(max_len)
                outname = record.id.split('.')[0] + ".ffn"
            rec_list.append(record)
        if genome_accession:
            header_accession = outname.split('.')[0]
    else:
        rec_list = seq_records
        header_accession = seq_records[0].name.split('.')[0]
    output_handle = open(outname, "w")
    for record in rec_list:
        for seq_feature in record.features:
            # skip pseudogenes 
            if 'pseudo' in seq_feature.qualifiers or 'pseudogene' in seq_feature.qualifiers:
                continue
            if seq_feature.type == "CDS":

                if locus_tag:
                    output_handle.write(">%s %s\n%s\n" % (
                            seq_feature.qualifiers["locus_tag"][0],
                            record.description,
                            seq_feature.extract(record.seq)))
                elif genome_accession:
                    output_handle.write(">%s %s\n%s\n" % (
                            header_accession,
                            record.description,
                            seq_feature.extract(record.seq)))

                else:
                    #print seq_feature
                    #try:
                    #    len(seq_feature.qualifiers['translation'])
                    #except:
                    #    print seq_feature
                    #    print "pseudogene?"
                    #    continue
                    #assert len(seq_feature.qualifiers['translation'])==1
                    # gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
                    try:
                        output_handle.write(">gi|%s|ref|%s| %s [%s]\n%s\n" % (
                                seq_feature.qualifiers["db_xref"][0].split(":")[1],
                                seq_feature.qualifiers["protein_id"][0],
                                seq_feature.qualifiers["product"][0],
                                record.description,
                                seq_feature.extract(record.seq)))
                    except:
                        try:

                            output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                                seq_feature.qualifiers["db_xref"][0].split(":")[1],
                                seq_feature.qualifiers["protein_id"][0],
                                #seq_feature.qualifiers["note"][0],
                                record.description,
                                seq_feature.extract(record.seq)))
                        except:
                            output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                                seq_feature.qualifiers["locus_tag"][0],
                                seq_feature.qualifiers["locus_tag"][0],
                                record.description,
                                seq_feature.extract(record.seq)))

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)
    parser.add_argument("-l", '--locus_tag', action='store_true', help="add locus_tag to header (optional)", default=False)
    parser.add_argument("-a", '--genome_accession', action='store_true', help="replace header of all CDS by the genome accession (optional)", default=False)

    args = parser.parse_args()



    input_handle = open(args.input_gbk, "rU")
    seq_records = list(SeqIO.parse(input_handle, "genbank"))

    if not args.outname:
        #print dir(seq_records[0])
        try:
            outname = seq_records[0].annotations['accessions'][0].split(".")[0]+".ffn"
        except:
            print(seq_records[0])
            import sys
            sys.exit()
        #print 'outname', outname
    else:
        outname = args.outname

    gbk2ffn(seq_records, outname, args.locus_tag, args.genome_accession)
