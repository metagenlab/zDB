#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def gbk2faa(seq_records,
            pformat=False,
            lformat=False,
            gformat=False,
            remove_redundancy=False,
            get_translation=False,
            outname=False,
            output_handle=False):

    import re, sys
    from Bio.SeqRecord import SeqRecord

    all_locus_ids = []
    all_prot_ids = []

    # hangeling of input str or input seqrecord
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
        outname = record_list[longest_record].id.split('.')[0] + ".faa"

    if not output_handle:
        output_handle = open(outname, "w")

    if pformat and lformat:
        raise('You have to chose either protein id or locus tag as header, not both!')


    '''
    for record in record_list:

        #input_handle = open(record, "rU")
        #seq_records = list(SeqIO.parse(input_handle, "genbank"))

    '''
    for record in record_list:
        description = record.description
        description = re.sub(", complete genome\.", "", description)
        description = re.sub(", complete sequence\.", "", description)
        description = re.sub("strain ", "", description)
        description = re.sub("str\. ", "", description)
        description = re.sub(" complete genome sequence\.", "", description)
        description = re.sub(" complete genome\.", "", description)
        description = re.sub(" chromosome", "", description)
        description = re.sub(" DNA", "S.", description)

        count_cds=1
        for seq_feature in record.features:
            if 'pseudo' in seq_feature.qualifiers or 'pseudogene' in seq_feature.qualifiers or not 'translation' in seq_feature.qualifiers:
                continue
            if seq_feature.type == "CDS":
                try:
                    if seq_feature.qualifiers["protein_id"][0] in all_prot_ids:
                        print ('%s (%s), Duplicated protein id: %s' % (record.id,
                                                                 record.description,
                                                                 seq_feature.qualifiers["protein_id"][0]))
                        if remove_redundancy:
                            print ('skipping...')
                            continue
                    else:
                        all_prot_ids.append(seq_feature.qualifiers["protein_id"][0])
                except KeyError:
                    # no protein ids
                    pass
                try:
                    if seq_feature.qualifiers["locus_tag"][0] in all_locus_ids:
                        print ('%s (%s), Duplicated locus id: %s, skipping' % (record.id,
                                                                 record.description,
                                                                 seq_feature.qualifiers["locus_tag"][0]))
                        # if remove redundancy, skip writing
                        if remove_redundancy:
                            print ('skipping...')
                            continue
                    else:
                        all_locus_ids.append(seq_feature.qualifiers["locus_tag"][0])
                except KeyError:
                    # no protein ids
                    pass
                # check presence of a protein sequence
                try:
                    len(seq_feature.qualifiers['translation'])
                except:
                    print (seq_feature)
                    print (record.seq)
                    import sys
                    sys.exit()
                    #sys.stderr.write(seq_feature.location.start)
                    sys.stderr.write("%s (%s) - %s: no sequence, is it a pseudogene, or genbank without translated CDS?\n" % (record.id,
                                                                 record.description,
                                                                 seq_feature.qualifiers["protein_id"][0]))
                    if get_translation:
                        import re
                        seq = str(seq_feature.extract(record.seq).translate())
                        if seq[-1] == '*':
                             seq = seq[0:-1]
                        seq_feature.qualifiers['translation'] = [seq]
                        print (seq)
                    else:
                        continue



                #assert len(seq_feature.qualifiers['translation'])==1
                # gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]

                if pformat:
                    # output only protein_id and species name
                    try:
                        output_handle.write(">%s %s\n%s\n" % (
                                seq_feature.qualifiers["protein_id"][0],
                                                description,
                                                seq_feature.qualifiers['translation'][0]))
                    except:
                        try:

                            output_handle.write(">%s %s\n%s\n" % (
                                seq_feature.qualifiers["protein_id"][0],
                                                description,
                                                seq_feature.qualifiers['translation'][0]))
                        except:
                            #print seq_feature
                            pass
                elif lformat:
                    # output only locus tag and species name
                    try:
                        output_handle.write(">%s %s\n%s\n" % (
                                seq_feature.qualifiers["locus_tag"][0],
                                                description,
                                                seq_feature.qualifiers['translation'][0]))
                    except:
                        try:

                            output_handle.write(">%s %s\n%s\n" % (
                                seq_feature.qualifiers["protein_id"][0],
                                                description,
                                                seq_feature.qualifiers['translation'][0]))
                        except:
                            #print seq_feature
                            pass
                elif gformat:
                    try:
                        output_handle.write(">%s_%s\n%s\n" % (
                                             record.name,
                                             count_cds,
                                             seq_feature.qualifiers['translation'][0]))
                        count_cds+=1
                    except:
                        print ('error', feature)
                    
                else:

                    try:
                        output_handle.write(">gi|%s|ref|%s| %s [%s]\n%s\n" % (
                                seq_feature.qualifiers["db_xref"][0].split(":")[1],
                                seq_feature.qualifiers["protein_id"][0],
                                seq_feature.qualifiers["product"][0],
                                description,
                                seq_feature.qualifiers['translation'][0]))
                    except:
                        try:

                            output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                                seq_feature.qualifiers["db_xref"][0].split(":")[1],
                                seq_feature.qualifiers["protein_id"][0],
                                #seq_feature.qualifiers["note"][0],
                                description,
                                seq_feature.qualifiers['translation'][0]))
                        except:
                            try:
                                output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                                    seq_feature.qualifiers["locus_tag"][0],
                                    seq_feature.qualifiers["locus_tag"][0],
                                    #seq_feature.qualifiers["note"][0],
                                    description,
                                    seq_feature.qualifiers['translation'][0]))
                            except:
                                if 'gene' in seq_feature.qualifiers:
                                    output_handle.write(">%s %s %s \n%s\n" % (
                                    seq_feature.qualifiers["gene"][0],
                                    seq_feature.qualifiers["product"][0],
                                    description,
                                    seq_feature.qualifiers['translation'][0]))
                                else:
                                    output_handle.write(">%s %s\n%s\n" % (
                                    seq_feature.qualifiers["product"][0],
                                    description,
                                    seq_feature.qualifiers['translation'][0]))



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file", nargs='+')
    #parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)
    parser.add_argument("-f", '--lformat', action='store_true', help="format header: >locus description", default=False)
    parser.add_argument("-p", '--pformat', action='store_true', help="format header: >protein_id description", default=False)
    parser.add_argument("-g", '--gformat', action='store_true', help="format header: >genome_accession", default=False)
    parser.add_argument("-o", '--outname', help="outname (optinal)", default=False)
    parser.add_argument("-r", '--remove', action='store_true', help="remove redundancy (protein id or locus tags persent mor than once)", default=False)
    parser.add_argument("-t", '--translate', action='store_true', help="translate from DNA if translation not available (not for pseudo tagged features)", default=False)
    parser.add_argument("-s", '--seqio', action='store_true',
                        help="simple SeqIO conversion",
                        default=False)
    args = parser.parse_args()

    if args.seqio:
        record_list = []
        for record in args.input_gbk:
            record_list += [i for i in SeqIO.parse(open(record), 'genbank')]
        with open(args.outname, 'w') as f:
            SeqIO.write(record_list, f, 'fasta')
    elif args.lformat:
        gbk2faa(args.input_gbk, lformat=args.lformat, remove_redundancy=args.remove, get_translation=args.translate, outname=args.outname)
    elif args.pformat:
        gbk2faa(args.input_gbk, pformat=args.pformat, remove_redundancy=args.remove, get_translation=args.translate, outname=args.outname)
    elif args.gformat:
        gbk2faa(args.input_gbk, gformat=args.gformat, remove_redundancy=args.remove, get_translation=args.translate, outname=args.outname)
    else:
        gbk2faa(args.input_gbk, False, False, remove_redundancy=args.remove, get_translation=args.translate, outname=args.outname)
