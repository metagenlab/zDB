#! /usr/bin/env python


import sys
from Bio import GenBank
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from chlamdb.biosqldb.manipulate_biosqldb import load_db
from chlamdb.biosqldb.manipulate_biosqldb import query_yes_no
from chlamdb.biosqldb import manipulate_biosqldb
import MySQLdb


def import_gbk(gbk_name,
               server_name,
               db_name):

    records = [i for i in SeqIO.parse(gbk_name, 'genbank')]

    if len(records)>1:
        print ("genbank file contains more than one record, is it plasmid(s) and chromosome(s)? (if the assembly contain more than one contig, contactenate contigs into a single record before importing it the the sql database)")
    try:
        db_name.load(records)
    except MySQLdb.Error as ee:
        print (gbk_name, "already into db?", ee)
        pass
    server_name.adaptor.commit()


def add_indexes(biodb_name):
    
    server, db = manipulate_biosqldb.load_db(biodb_name)
    
    sql_index1 = 'create index ftgcga on feature_tables_genomes_cds(genome_accession)'    
    sql_index2 = 'create index ftgctx on feature_tables_genomes_cds(taxon_id)'
    sql_index3 = 'create index ftgrga on feature_tables_genomes_rrna(taxon_id)'
    sql_index4 = 'create index ftgrtx on feature_tables_genomes_rrna(genome_accession)'    
    sql_index5 = 'create index ftcain on feature_tables_cds_accessions(id_name)'
    sql_index6 = 'create index ftcait on feature_tables_cds_accessions(id_type)'
    
    server.adaptor.execute(sql_index1,)
    server.adaptor.execute(sql_index2,)
    server.adaptor.execute(sql_index3,)
    server.adaptor.execute(sql_index4,)
    server.adaptor.execute(sql_index5,)
    server.adaptor.execute(sql_index6,)
    
def create_cds_tables(one_gbk,
                      biodb_name):

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

    import re
    server, db = manipulate_biosqldb.load_db(biodb_name)

    # TODO add unique constrains to prevent inserting 2 times the same genome

    sql_cds = 'CREATE table IF NOT EXISTS feature_tables_genomes_cds (prot_primary_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' taxon_id INT,' \
           ' genome_accession VARCHAR(40),' \
           ' start INT,' \
           ' end INT,' \
           ' strand INT,' \
           ' gene varchar(20),' \
           ' product TEXT,' \
           ' translation TEXT)'
    

    
    sql_rrna = 'CREATE table IF NOT EXISTS feature_tables_genomes_rrna (rrna_primary_id INT AUTO_INCREMENT PRIMARY KEY,' \
               ' taxon_id INT,' \
               ' genome_accession VARCHAR(40),' \
               ' start INT,' \
               ' end INT,' \
               ' strand INT,' \
               ' product TEXT)'


    sql_synonyms = 'CREATE table IF NOT EXISTS feature_tables_cds_accessions (prot_primary_id INT,' \
                   ' id_type varchar(40),' \
                   ' id_name varchar(40))' \
                   #' FOREIGN KEY (prot_primary_id) REFERENCES feature_tables_genomes_cds(prot_primary_id))'
                   

    
    server.adaptor.execute(sql_cds,)
    server.adaptor.execute(sql_rrna,)
    server.adaptor.execute(sql_synonyms,)
    
    
    server.commit()
    

    with open(one_gbk, 'r') as f:
        records = SeqIO.parse(f, 'genbank')
        for record in records:
            accession = record.id.split('.')[0]
            
            sql = 'select taxon_id from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                  ' where t2.name="%s" and t1.accession="%s"' % (biodb_name,
                                                                 accession)

            taxon_id = server.adaptor.execute_and_fetchall(sql,)[0][0]
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

                    try:
                        protein_id = feature.qualifiers['protein_id'][0]
                    except:
                        protein_id = '.'
                    try:
                        translation = feature.qualifiers['translation'][0]
                    except:
                        translation = '-'

                    sql1 = 'INSERT INTO feature_tables_genomes_cds(taxon_id, genome_accession, start, end, strand, gene, product, translation) ' \
                           'values(%s, "%s", %s, %s, %s, "%s", "%s", "%s")' % (taxon_id,
                                                                               accession,
                                                                               start,
                                                                               end,
                                                                               strand,
                                                                               gene,
                                                                               product,
                                                                               translation)
                    server.adaptor.execute(sql1,)
                    server.commit()

                    cds_id = server.adaptor.cursor.lastrowid
                    sql2 = 'INSERT into feature_tables_cds_accessions(prot_primary_id, id_type, id_name) values (' \
                           ' %s, "%s", "%s")'
                    server.adaptor.execute(sql2 % (cds_id, 'protein_id', protein_id),)
                    server.adaptor.execute(sql2 % (cds_id, 'locus_tag', locus_tag),)

                    if old_locus_tag:
                        server.adaptor.execute(sql2 % (cds_id, 'old_locus_tag', old_locus_tag),)

                    server.commit()

                elif feature.type == 'rRNA':
                    start = re.sub('>|<','', str(feature.location.start))
                    end = re.sub('>|<','', str(feature.location.end))
                    strand = feature.strand
                    try:
                        product = feature.qualifiers['product'][0]
                    except:
                        product = '-'
                    sql = 'insert into feature_tables_genomes_rrna (taxon_id, genome_accession, start, end, strand, product) values (' \
                          ' %s, "%s", %s, %s, %s, "%s")' % (taxon_id, accession, start, end, strand, product)

                    server.adaptor.execute(sql,)
                    server.commit()
                else:
                    continue

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", '--gbk', type=str, help="gbk file", nargs='+')
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-r", '--remove_db', type=str, help="remove db")

    args = parser.parse_args()

    server, db = load_db(args.db_name)

    if args.gbk:
        sys.stdout.write("Importing Gbk files & Creating CDS tables...\n")
        for gbk in args.gbk:
            print(gbk)
            import_gbk(gbk, server, db)
            create_cds_tables(gbk, args.db_name)
        sys.stdout.write("Update config table...\n")
        manipulate_biosqldb.update_config_table(args.db_name, "gbk_files")
        try:
            add_indexes(args.db_name)
        except:
            print("Tables already indexed?")

    if args.remove_db:
        if query_yes_no("Remove databse %s?" % args.remove_db):
            del server[args.remove_db]
            server.commit()
