#! /usr/bin/env python


import sys
from Bio import GenBank
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from manipulate_biosqldb import load_db
from manipulate_biosqldb import query_yes_no

def import_gbk(gbk_name, server_name, db_name, sqlite=False):

    from Bio import SeqIO
    #parser = GenBank.FeatureParser()
    #iterator = GenBank.Iterator(open(gbk_name), parser)
    records = [i for i in SeqIO.parse(gbk_name, 'genbank')]
    #iterator = GenBank.Iterator(open(gbk_name), parser)
    if len(records)>1:
        print "genbank file contains more than one record, is it plasmid(s) and chromosome(s)? (if the assembly contain more than one contig, contactenate contigs into a single record before importing it the the sql database)"
        #sys.exit()
    #With this iterator, the loading of the database is another one-liner:
    db_name.load(records)
    try:
        db_name.load(records)
    except:
        db_name.load(records)
        print gbk_name, "already into db?"
        pass
    server_name.adaptor.commit()




def import_swiss(swiss_dat_file, server_name, db):

    db.load(SeqIO.parse(swiss_dat_file, "swiss"))
    server_name.adaptor.commit()


def create_cds_tables(one_gbk, biodb_name):

    '''

    Create one protein CDS and one rRNA table for each genbank record.
    The genbak should already be loaded into the biosql database.
    Columns of the protein encoding cds table:

    - taxon_id index
    - accession index
    - locus_tag index
    - protein_id index
    - start
    - stop
    - strand
    - gene
    - product
    - translation

    :param one_gbk:
    :return:
    '''
    import manipulate_biosqldb
    import re
    server, db = manipulate_biosqldb.load_db(biodb_name)

    sql_cds = 'CREATE table IF NOT EXISTS feature_tables.genomes_cds_%s (prot_primary_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' taxon_id INT,' \
           ' genome_accession VARCHAR(40),' \
           ' start INT,' \
           ' end INT,' \
           ' strand INT,' \
           ' gene varchar(20),' \
           ' product TEXT,' \
           ' translation TEXT, ' \
           ' INDEX taxon_id(taxon_id), ' \
           ' INDEX genome_accession(genome_accession))' % biodb_name

    sql_rrna = 'CREATE table IF NOT EXISTS feature_tables.genomes_rrna_%s (rrna_primary_id INT AUTO_INCREMENT PRIMARY KEY,' \
               ' taxon_id INT,' \
               ' genome_accession VARCHAR(40),' \
               ' start INT,' \
               ' end INT,' \
               ' strand INT,' \
               ' product TEXT,' \
               ' INDEX taxon_id(taxon_id), ' \
               ' INDEX genome_accession(genome_accession))' % biodb_name

    sql_synonyms = 'CREATE table IF NOT EXISTS feature_tables.cds_accessions_%s (prot_primary_id INT,' \
                   ' id_type varchar(40),' \
                   ' id_name varchar(40),' \
                   ' INDEX id_name(id_name),' \
                   ' INDEX id_type(id_type),' \
                   ' FOREIGN KEY (prot_primary_id) REFERENCES genomes_cds_%s(prot_primary_id))' % (biodb_name, biodb_name)

    print sql_cds
    server.adaptor.execute(sql_cds,)
    print sql_rrna
    server.adaptor.execute(sql_rrna,)
    print sql_synonyms
    server.adaptor.execute(sql_synonyms,)
    server.commit()


    with open(one_gbk, 'r') as f:
        records = SeqIO.parse(f, 'genbank')
        for record in records:
            accession = record.id.split('.')[0] # fill_color=violet
            sql = 'select taxon_id from bioentry as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                  ' where t2.name="%s" and t1.accession="%s"' % (biodb_name,
                                                                 accession)
            taxon_id = server.adaptor.execute_and_fetchall(sql,)[0][0]
            print 'taxon id', taxon_id

            for feature in record.features:
                if feature.type == 'CDS' and not 'pseudo' in feature.qualifiers:
                    start = re.sub('>|<','', str(feature.location.start))
                    end = re.sub('>|<','', str(feature.location.end))
                    strand = feature.strand
                    try:
                        gene = feature.qualifiers['gene'][0]
                    except:
                        gene = '-'
                    try:
                        product = feature.qualifiers['product'][0]
                    except:
                        product = '-'
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    try:
                        old_locus_tag = feature.qualifiers['old_locus_tag'][0]
                    except:
                        old_locus_tag = False

                    protein_id = feature.qualifiers['protein_id'][0]
                    translation = feature.qualifiers['translation'][0]


                    sql1 = 'INSERT INTO feature_tables.genomes_cds_%s(taxon_id, genome_accession, start, end, strand, gene, product, translation) ' \
                           'values(%s, "%s", %s, %s, %s, "%s", "%s", "%s")' % (biodb_name,
                                                                               taxon_id,
                                                                               accession,
                                                                               start,
                                                                               end,
                                                                               strand,
                                                                               gene,
                                                                               product,
                                                                               translation)
                    server.adaptor.execute(sql1,)
                    server.commit()

                    sql = 'SELECT LAST_INSERT_ID();'

                    cds_id = server.adaptor.execute_and_fetchall(sql, )[0][0]
                    sql2 = 'INSERT into feature_tables.cds_accessions_%s(prot_primary_id, id_type, id_name) values (' \
                           ' %s, "%s", "%s")'

                    server.adaptor.execute(sql2 % (biodb_name, cds_id, 'protein_id', protein_id),)
                    server.adaptor.execute(sql2 % (biodb_name, cds_id, 'locus_tag', locus_tag),)

                    if old_locus_tag:
                        server.adaptor.execute(sql2 % (biodb_name, cds_id, 'old_locus_tag', old_locus_tag),)

                    server.commit()

                elif feature.type == 'rRNA':
                    start = re.sub('>|<','', str(feature.location.start))
                    end = re.sub('>|<','', str(feature.location.end))
                    strand = feature.strand
                    product = feature.qualifiers['product'][0]

                    sql = 'insert into feature_tables.genomes_rrna_%s (taxon_id, genome_accession, start, end, strand, product) values (' \
                          ' %s, "%s", %s, %s, %s, "%s")' % (biodb_name, taxon_id, accession, start, end, strand, product)

                    server.adaptor.execute(sql,)
                    server.commit()
                else:
                    continue

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", '--gbk', type=str, help="gbk file", nargs='+')
    parser.add_argument("-u", '--uniprot', type=str, help="uniprot dat file", nargs='+')
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-r", '--remove_db', type=str, help="remove db")
    parser.add_argument("-s", '--sqlite_db', type=str, help="use sqlite3 database", default=False)

    
    args = parser.parse_args()


    server, db = load_db(args.db_name, args.sqlite_db)

    if args.gbk:
        for gbk in args.gbk:
            print gbk
            try:
                import_gbk(gbk, server, db)
            except:
                print 'problem importing %s' % gbk
                #import_gbk(gbk, server, db)
            #create_cds_tables(gbk, args.db_name)
    
    if args.remove_db:
        if query_yes_no("Remove databse %s?" % args.remove_db):
            del server[args.remove_db]
            server.commit()

    if args.uniprot:
        for swiss in args.uniprot:
            import_swiss(swiss, server, db)
