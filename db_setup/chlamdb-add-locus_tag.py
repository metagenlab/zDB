#!/usr/bin/env python


def add_missing_locus__tag(one_gbk):

    from Bio import SeqIO
    records = [i for i in SeqIO.parse(one_gbk, 'genbank')]

    for i, record in enumerate(records):
        protein2count = {}
        for n, feature in enumerate(record.features):
            if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers and 'pseudogene' not in feature.qualifiers:
                if not 'locus_tag' in feature.qualifiers:
                    print("Adding locus tag:\t%s" % (feature.qualifiers["protein_id"][0]))
                    protein_id = feature.qualifiers["protein_id"][0].split(".")[0]
                    if protein_id not in protein2count:
                        protein2count[protein_id] = 1
                        locus_tag = protein_id
                    else:
                        protein2count[protein_id] += 1
                        locus_tag = "%s_%s" % (protein_id, protein2count[protein_id])
                    records[i].features[n].qualifiers["locus_tag"] = [locus_tag]
    SeqIO.write(records, "%s_locus.gbk" % one_gbk.split(".")[0], "genbank")


if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", '--gbk_file', type=str, help="GBK file")

    args = parser.parse_args()

    add_missing_locus__tag(args.gbk_file)
