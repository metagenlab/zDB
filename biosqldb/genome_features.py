#!/usr/bin/python

import manipulate_biosqldb
def get_n_cds(server, biodb_name, bioentry_id):
    sql = 'select count(t1.value) from biosqldb.seqfeature_qualifier_value as t1' \
          ' inner join seqfeature as t2 on t1.seqfeature_id = t2.seqfeature_id'\
        ' inner join term as t3 on t1.term_id = t3.term_id and t3.name = "translation"' \
        ' inner join bioentry as t6 on t2.bioentry_id = t6.bioentry_id and t6.accession = "%s"' \
        ' inner join biodatabase as t7 on t6.biodatabase_id = t7.biodatabase_id and t7.name = "%s"' % (bioentry_id, biodb_name)
    result = server.adaptor.execute_and_fetchall(sql, )

    return int(result[0][0])


def get_n_contigs(db, accession):
    n_contigs = 1
    genome = db.lookup(accession=accession)
    for feature in genome.features:
        if feature.type=="assembly_gap":
            n_contigs+=1
    return n_contigs


if __name__ == '__main__':
    server, db = manipulate_biosqldb.load_db("chlamydia_03_15")

    print get_n_cds(server, "chlamydia_03_15", "CP001848")
    accessions = ['ACZE01000000', 'APJW01000000', 'AYKJ01000001', 'CCEJ010000000', 'CCSC01000000', 'JSAM01000000', 'JSAN01000000', 'JSDQ01000000', 'JRXI01000000', 'NC_004552', 'AE001273', 'AE015925', 'AE015926', 'AYKJ01000000', 'NC_007899', 'NC_007900', 'NC_002620', 'AE002162', 'CP001928', 'CP001929', 'CP002549', 'CP002550', 'CP002608', 'CP006571', 'CP006572', 'NC_002179', 'ElaC', 'ElaCp', 'FR872580', 'CP000975', 'KNic', 'KNicp', 'BASK00000000', 'NHA', 'BAWW00000000', 'NZ_BBPT00000000', 'NC_005861', 'PnaD', 'NZ_BASL01000000', 'CP001848', 'BX119912', 'Rhab', 'Rhabp', 'NC_015713', 'NC_015710']
    n_contigs = []

    for i in accessions:
        print i
        n_contigs.append(get_n_contigs(db, i))
    print n_contigs