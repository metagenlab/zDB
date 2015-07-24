#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# create html table from blast tabulated file
# headers: accession	size	gi	n proteins	n contigs	gc 	description
# add 4 columns with links of the form /assets/chlamdb/ffn/ for gbk/faa/ffm/fna
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2015
# ---------------------------------------------------------------------------


def interpro2biosql(server,
                    seqfeature_id2locus_tag,
                    locus_tag2genome_taxon_id,
                    protein_id2genome_taxon_id,
                    locus_tag2seqfeature_id,
                    protein_id2seqfeature_id,
                    seqfeature_id2organism,
                    db_name,
                    *input_files):
    import re
    '''
    1. Protein Accession (e.g. P51587)
    2. Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    3.. equence Length (e.g. 3418)
    4. Analysis (e.g. Pfam / PRINTS / Gene3D)
    5. Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
    6. Signature Description (e.g. BRCA2 repeat profile)
    7. Start location
    8. Stop location
    9. Score - is the e-value of the match reported by member database method (e.g. 3.1E-52)
    10. Status - is the status of the match (T: true)
    11. Date - is the date of the run
    12. (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprscan option is switched on)
    13. (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprscan option is switched on)
    14. (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
    15. (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
    :param input_file:
    :return:
    '''

    sql = 'CREATE TABLE interpro_%s (accession VARCHAR(100),' \
          ' locus_tag VARCHAR(200), ' \
          ' organism VARCHAR(200),  ' \
          ' taxon_id INT,' \
          ' sequence_length INT, ' \
          ' analysis VARCHAR(100) NOT NULL, ' \
          ' signature_accession VARCHAR(100), ' \
          ' signature_description VARCHAR(1000), ' \
          ' start INT, ' \
          ' stop INT, ' \
          ' score VARCHAR(10) NOT NULL, ' \
          ' interpro_accession VARCHAR(1000) NOT NULL, ' \
          ' interpro_description VARCHAR(10000),' \
          ' GO_terms varchar(10000),' \
          ' pathways varchar(10000))' % db_name
    try:
        server.adaptor.execute(sql)
    except:
        pass
    for one_interpro_file in input_files:
        print one_interpro_file
        from pandas import read_csv
        with open(one_interpro_file, 'r') as f:
            tsvin = read_csv(f, sep='\t', error_bad_lines=False, names=["accession",
                                                                        "MD5",
                                                                        "sequence_length",
                                                                        "analysis",
                                                                        "signature_accession",
                                                                        "signature_description",
                                                                        "start",
                                                                        "stop",
                                                                        "score",
                                                                        "status",
                                                                        "date",
                                                                        "interpro_accession",
                                                                        "interpro_description",
                                                                        "GO_terms",
                                                                        "pathways"])
            tsvin = tsvin.fillna(0)

            for i in range(len(tsvin['accession'])):

                data= list(tsvin.loc[i,:])
                for index, item in enumerate(data):
                    if type(item) == str:
                        data[index] = item.replace('\"','')

                accession = data[0]

                sequence_length = data[2]
                analysis = data[3]
                signature_accession = data[4]
                signature_description = data[5]
                start = data[6]
                stop = data[7]
                score = data[8]


                interpro_accession = data[11]
                interpro_description = data[12]
                GO_terms = data[13]
                pathways = data[14]

                try:
                    taxon_id = protein_id2genome_taxon_id[accession]
                    seqfeature_id = protein_id2seqfeature_id[accession]
                except KeyError:
                    taxon_id = locus_tag2genome_taxon_id[accession]
                    seqfeature_id = locus_tag2seqfeature_id[accession]
                organism = seqfeature_id2organism[str(seqfeature_id)]
                locus_tag = seqfeature_id2locus_tag[str(seqfeature_id)]
                sql = 'INSERT INTO interpro_%s(accession, locus_tag, organism, taxon_id,' \
                      ' sequence_length, analysis, signature_accession, signature_description, start, ' \
                      ' stop, score, interpro_accession, interpro_description, GO_terms, pathways) ' \
                      ' values ("%s", "%s", "%s", %s, %s, "%s", "%s", "%s", %s, %s, "%s", "%s", "%s", "%s", "%s");' % (db_name,
                                                                                                     accession,
                                                                                                     locus_tag,
                                                                                                     organism,
                                                                                                     taxon_id,
                                                                                               sequence_length,
                                                                                               analysis,
                                                                                               signature_accession,
                                                                                               signature_description,
                                                                                               int(start),
                                                                                               int(stop),
                                                                                               str(score),
                                                                                               interpro_accession,
                                                                                               interpro_description,
                                                                                               GO_terms,
                                                                                               pathways)
                try:
                    server.adaptor.execute(sql)
                    server.adaptor.commit()
                except:
                    print sql
                    print data
                    import sys
                    sys.exit()
if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_interpro', type=str, help="input interpro xml file", nargs='+')

    args = parser.parse_args()

    biodb = 'chlamydia_03_15'

    server, db = manipulate_biosqldb.load_db(biodb)
    print "creating locus_tag2seqfeature_id"
    locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)
    print "creating protein_id2seqfeature_id"
    protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, biodb)

    print "getting seqfeature_id2organism"
    seqfeature_id2organism = manipulate_biosqldb.seqfeature_id2organism_dico(server, biodb)

    print "creating locus_tag2taxon_id dictionnary..."
    locus_tag2genome_taxon_id = manipulate_biosqldb.locus_tag2genome_taxon_id(server, biodb)
    print "creating protein_id2taxon_id dictionnary..."
    protein_id2genome_taxon_id = manipulate_biosqldb.protein_id2genome_taxon_id(server, biodb)
    print "getting seqfeature_id2locus_tag"
    seqfeature_id2locus_tag = manipulate_biosqldb.seqfeature_id2locus_tag_dico(server, biodb)


    interpro2biosql(server,seqfeature_id2locus_tag, locus_tag2genome_taxon_id, protein_id2genome_taxon_id, locus_tag2seqfeature_id, protein_id2seqfeature_id, seqfeature_id2organism, biodb, *args.input_interpro)


