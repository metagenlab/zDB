#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# Load blast results into biosql databas
# Tablulated blast files: qgi qacc sgi sacc sscinames sskingdoms staxids evalue nident pident positive gaps length qstart qend qcovs sstart send sstrand stitle
# Genbank files of query proteins must be leaded in the biosql database already (used to fetch locus tag, taxon_id, organism name, ...)
# Create 3 new tables (if they do not exit): - blastnr_taxonomy with txaon_id and corresponding full taxonomic path
#                                            todo: use taxonomy implemented in biosqldb shema instead: currently fetch taxonomic path from ncbi and add it to blast_nr_taxonomy
#                                           - blast_nr_`biodb_name`: contain data foreach hit
#                                           - blast_nr_taxonomy_`biodb_name`: contain taxonomical data oreach hit (as one hit can now include multiple taxons/MULTISPECIES hits)
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: april 2015
# ---------------------------------------------------------------------------

'''

    #### general data ####

    hit_number
    taxon_id
    locus_tag
    organism
    orthogroup

    #### NCBI taxonomical ranks ####

    query_gi
    subject_gi
    subject_scientific_name
    subject_taxid
    evalue
    n_identical
    percent_identity
    positive
    gaps
    length
    query_start
    query_end
    query_cov
    subject_start
    subject_end]
    subject_strand
    subject_title

    no_rank
    superkingdom
    kingdom
    subkingdom
    superphylum
    phylum
    subphylum
    superclass
    class
    subclass
    superorder
    order
    suborder
    superfamily
    family
    subfamily
    genus
    subgenus
    species
    species_subgroup
    species_group
    subspecies

    tribe
    infraorder
    subtribe
    forma
    infraclass
    varietas
    parvorder
'''



def blastnr2biosql(server,
                    seqfeature_id2locus_tag,
                    locus_tag2genome_taxon_id,
                    protein_id2genome_taxon_id,
                    locus_tag2seqfeature_id,
                    protein_id2seqfeature_id,
                    seqfeature_id2organism,
                    db_name, create_tables=False, *input_blast_files):

    import sequence_id2scientific_classification
    import time

    sql_blast = 'CREATE TABLE IF NOT EXISTS blastnr_%s (nr_hit_id INT AUTO_INCREMENT PRIMARY KEY, ' \
                ' query_accession varchar(200),' \
                ' hit_number int,' \
                ' taxon_id int,' \
                ' locus_tag VARCHAR(100), ' \
                ' organism VARCHAR(100), ' \
                ' query_gi int,' \
                ' subject_gi int,' \
                ' subject_accession varchar(200),' \
                ' subject_kingdom varchar(200),' \
                ' subject_scientific_name TEXT(2000000), ' \
                ' subject_taxid TEXT(2000000),' \
                ' evalue float,' \
                ' n_identical int,' \
                ' percent_identity float,' \
                ' positive int,' \
                ' gaps int,' \
                ' length int,' \
                ' query_start int,' \
                ' query_end int,' \
                ' query_cov int,' \
                ' subject_start int,' \
                ' subject_end int,' \
                ' subject_strand VARCHAR(5), ' \
                ' subject_title VARCHAR(2000))' % db_name

    sql_taxonomy1 = 'CREATE TABLE IF NOT EXISTS blastnr_taxonomy (taxon_id int unique, ' \
                    ' no_rank VARCHAR(200) default "-", ' \
                    ' superkingdom VARCHAR(200) default "-", ' \
                    ' kingdom VARCHAR(200) default "-", ' \
                    ' subkingdom VARCHAR(200) default "-", ' \
                    ' superphylum VARCHAR(200) default "-", ' \
                    ' phylum VARCHAR(200) default "-", ' \
                    ' subphylum VARCHAR(200) default "-", ' \
                    ' superclass VARCHAR(200) default "-", ' \
                    ' class VARCHAR(200) default "-", ' \
                    ' subclass VARCHAR(200) default "-", ' \
                    ' superorder VARCHAR(200) default "-", ' \
                    ' `order` VARCHAR(200) default "-", ' \
                    ' suborder VARCHAR(200) default "-", ' \
                    ' superfamily VARCHAR(200) default "-", ' \
                    ' family VARCHAR(200) default "-", ' \
                    ' subfamily VARCHAR(200) default "-", ' \
                    ' genus VARCHAR(200) default "-", ' \
                    ' subgenus VARCHAR(200) default "-", ' \
                    ' species VARCHAR(200) default "-", ' \
                    ' species_subgroup VARCHAR(200) default "-", ' \
                    ' species_group VARCHAR(200) default "-", ' \
                    ' subspecies VARCHAR(200) default "-", ' \
                    ' tribe VARCHAR(200) default "-", ' \
                    ' infraorder VARCHAR(200) default "-", ' \
                    ' subtribe VARCHAR(200) default "-", ' \
                    ' forma VARCHAR(200) default "-", ' \
                    ' infraclass VARCHAR(200) default "-", ' \
                    ' varietas VARCHAR(200) default "-", ' \
                    ' parvorder VARCHAR(200) default "-")'


    sql_taxonomy2 = 'CREATE TABLE IF NOT EXISTS blastnr_taxonomy_%s (nr_hit_id INT, ' \
                    ' hit_number int(10),' \
                    ' taxon_id int)' % (db_name)

    sql_taxonomy3 = 'ALTER TABLE blastnr_taxonomy_%s ADD CONSTRAINT fk_blast_nr_hit_id FOREIGN KEY (nr_hit_id) REFERENCES blastnr_%s(nr_hit_id);' % (db_name, db_name)

    #try:
    #print sql_blast
    if create_tables:
        server.adaptor.execute(sql_blast)
        #print sql_taxonomy
        server.adaptor.execute(sql_taxonomy1)
        server.adaptor.execute(sql_taxonomy2)
        server.adaptor.execute(sql_taxonomy3)
        #server.adaptor.commit()
        #except:
        #    print 'blastnr tables already exist!'
    print 'all blast files'
    print input_blast_files
    for one_blast_file in input_blast_files:
                print 'blast file %s' % one_blast_file
                with open(one_blast_file, 'r') as f:
                    input_file = [i.rstrip().split('\t') for i in f]
                    hit_n = 0
                    sql = 'SELECT taxon_id from blastnr_taxonomy'
                    all_taxon_ids = []
                    database_taxon_ids = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
                    print "Number of taxons into database: ", len(database_taxon_ids)
                    for n, line in enumerate(input_file):
                        temp_subject_taxids = line[6].split(';')
                        # store taxon_ids
                        for i in temp_subject_taxids:
                            if i not in database_taxon_ids and i not in all_taxon_ids:
                                all_taxon_ids.append(i)
                    print 'Fetching ncbi taxonomy... for %s taxons' % str(len(all_taxon_ids))
                    taxid2classification = sequence_id2scientific_classification.taxon_id2scientific_classification(all_taxon_ids)

                    print 'Updating blastnr_taxonomy table...'
                    for taxon_id in all_taxon_ids:
                            if taxon_id == 'N/A':
                                continue
                            import re
                            import MySQLdb
                            sql = 'INSERT INTO blastnr_taxonomy(taxon_id) values (%s)' % (taxon_id)
                            try:

                                server.adaptor.execute(sql)
                                server.adaptor.commit()
                            except MySQLdb.IntegrityError:
                                print 'Taxon %s already in database' % str(taxon_id)
                                continue
                            for rank in taxid2classification[taxon_id].keys():

                                sql_id = re.sub(' ', '_', rank)

                                value = taxid2classification[taxon_id][rank]
                                sql = 'UPDATE blastnr_taxonomy SET `%s`="%s" where taxon_id=%s' % (sql_id, value, taxon_id)
                                #print sql

                                server.adaptor.execute(sql)
                                server.adaptor.commit()

                            '''

                            sql_taxonomy_insert =   'INSERT INTO blastnr_taxonomy(taxon_id, ' \
                                                    'no_rank, ' \
                                                    'superkingdom, ' \
                                                    'kingdom, ' \
                                                    'subkingdom, ' \
                                                    'superphylum, ' \
                                                    'phylum, ' \
                                                    'subphylum, ' \
                                                    'superclass, ' \
                                                    'class, ' \
                                                    'subclass, ' \
                                                    'superorder, ' \
                                                    'taxorder, ' \
                                                    'suborder, ' \
                                                    'superfamily, ' \
                                                    'family, ' \
                                                    'subfamily, ' \
                                                    'genus, ' \
                                                    'subgenus, ' \
                                                    'species, ' \
                                                    'species_subgroup, ' \
                                                    'species_group, ' \
                                                    'subspecies, ' \
                                                    'tribe, ' \
                                                    'infraorder, ' \
                                                    'subtribe, ' \
                                                    'forma, ' \
                                                    'infraclass, ' \
                                                    'varietas, ' \
                                                    'parvorder)' \
                                                    ' values (%s, "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s",' \
                                                    '"%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s",' \
                                                    ' "%s", "%s");' % (taxon_id,
                                                       taxid2classification[taxon_id]['no rank'],
                                                        taxid2classification[taxon_id]['superkingdom'],
                                                        taxid2classification[taxon_id]['kingdom'],
                                                        taxid2classification[taxon_id]['subkingdom'],
                                                        taxid2classification[taxon_id]['superphylum'],
                                                        taxid2classification[taxon_id]['phylum'],
                                                        taxid2classification[taxon_id]['subphylum'],
                                                        taxid2classification[taxon_id]['superclass'],
                                                        taxid2classification[taxon_id]['class'],
                                                        taxid2classification[taxon_id]['subclass'],
                                                        taxid2classification[taxon_id]['superorder'],
                                                        taxid2classification[taxon_id]['order'],
                                                        taxid2classification[taxon_id]['suborder'],
                                                        taxid2classification[taxon_id]['superfamily'],
                                                        taxid2classification[taxon_id]['family'],
                                                        taxid2classification[taxon_id]['subfamily'],
                                                        taxid2classification[taxon_id]['genus'],
                                                        taxid2classification[taxon_id]['subgenus'],
                                                        taxid2classification[taxon_id]['species'],
                                                        taxid2classification[taxon_id]['species_subgroup'],
                                                        taxid2classification[taxon_id]['species_group'],
                                                        taxid2classification[taxon_id]['subspecies'],
                                                        taxid2classification[taxon_id]['tribe'],
                                                        taxid2classification[taxon_id]['infraorder'],
                                                        taxid2classification[taxon_id]['subtribe'],
                                                        taxid2classification[taxon_id]['forma'],
                                                        taxid2classification[taxon_id]['infraclass'],
                                                        taxid2classification[taxon_id]['varietas'],
                                                        taxid2classification[taxon_id]['parvorder']
                                                        )

                            '''

                    print 'loading blast results into database...'
                    for n, line in enumerate(input_file):
                        # qgi qacc sgi sacc sscinames sskingdoms staxids evalue nident pident positive gaps length qstart qend qcovs sstart send sstrand stitle

                        if n%1000 == 0:
                            print time.ctime() + ': ' + round((float(n)/len(input_file))*100,2), '%...'
                        query_accession = line[1].split("|")[3]

                        if n == 0:
                            hit_n+=1
                            pass
                        elif line[1] != input_file[n-1][1]:
                            hit_n = 1
                        else:
                            # handleling of multiple hsp
                            if line[3] != input_file[n-1][3]:
                                hit_n+=1
                            pass

                        subject_accession = line[3]
                        query_gi = int(line[0])
                        subject_gi = int(line[2])

                        temp_subject_scientific_names = line[4].split(';')

                        if len(temp_subject_scientific_names) == 1:
                            subject_scientific_names = temp_subject_scientific_names[0]
                        else:
                            subject_scientific_names = ''
                            for i in range(0, len(temp_subject_scientific_names)-1):
                                subject_scientific_names += '%s, ' % temp_subject_scientific_names[i]
                            subject_scientific_names += temp_subject_scientific_names[-1]

                        subject_kingdom = line[5]
                        temp_subject_taxids = line[6].split(';')


                        if len(temp_subject_taxids) == 1:
                            subject_taxids = temp_subject_taxids[0]
                        else:
                            subject_taxids = ''
                            for i in range(0, len(temp_subject_taxids)-1):
                                subject_taxids += '%s, ' % temp_subject_taxids[i]
                            subject_taxids += temp_subject_taxids[-1]



                        evalue = float(line[7])
                        n_identical = int(line[8])
                        percent_identity = float(line[9])
                        positive = int(line[10])
                        gaps = int(line[11])
                        length = int(line[12])
                        query_start = int(line[13])
                        query_end = int(line[14])
                        query_cov = int(line[15])
                        subject_start = int(line[16])
                        subject_end = int(line[17])
                        subject_strand = line[18]
                        subject_title = line[19]



                        try:
                            taxon_id = protein_id2genome_taxon_id[query_accession]
                            seqfeature_id = protein_id2seqfeature_id[query_accession]
                        except KeyError:
                            taxon_id = locus_tag2genome_taxon_id[query_accession]
                            seqfeature_id = locus_tag2seqfeature_id[query_accession]
                        organism = seqfeature_id2organism[str(seqfeature_id)]
                        locus_tag = seqfeature_id2locus_tag[str(seqfeature_id)]

                        sql_blast = 'INSERT INTO blastnr_%s(hit_number, ' \
                                    'query_accession, ' \
                                'taxon_id, ' \
                                'locus_tag, ' \
                                'organism, ' \
                                'query_gi, ' \
                                'subject_gi, ' \
                                'subject_accession, ' \
                                'subject_kingdom, ' \
                                'subject_scientific_name, ' \
                                'subject_taxid, ' \
                                'evalue, ' \
                                'n_identical, ' \
                                'percent_identity, ' \
                                'positive, ' \
                                'gaps, ' \
                                'length, ' \
                                'query_start, ' \
                                'query_end, ' \
                                'query_cov, ' \
                                'subject_start, ' \
                                'subject_end, ' \
                                'subject_strand, ' \
                                'subject_title)' \
                                ' values (%s,"%s", %s,"%s","%s",%s,%s,"%s","%s","%s","%s",%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,"%s","%s")' % (db_name,
                                                                                                                        hit_n,
                                                                                                                        query_accession,
                                                                                                                        taxon_id,
                                                                                                                        locus_tag,
                                                                                                                        organism,
                                                                                                                        query_gi,
                                                                                                                        subject_gi,
                                                                                                                        subject_accession,
                                                                                                                        subject_kingdom,
                                                                                                                        subject_scientific_names,
                                                                                                                        subject_taxids,
                                                                                                                        evalue,
                                                                                                                        n_identical,
                                                                                                                        percent_identity,
                                                                                                                        positive,
                                                                                                                        gaps,
                                                                                                                        length,
                                                                                                                        query_start,
                                                                                                                        query_end,
                                                                                                                        query_cov,
                                                                                                                        subject_start,
                                                                                                                        subject_end,
                                                                                                                        subject_strand,
                                                                                                                        subject_title
                                                                                                                        )


                        try:
                            server.adaptor.execute(sql_blast)
                            server.adaptor.commit()
                        except:
                            print 'problem with:'
                            print sql_blast
                            print line
                            import sys
                            sys.exit()
                        sql1 = 'select nr_hit_id from blastnr_%s where (locus_tag="%s" and hit_number=%s)' % (db_name,
                                                                                                              locus_tag,
                                                                                                              hit_n)
                        blast_id = server.adaptor.execute_and_fetchall(sql1,)[0][0]

                        # new table for the multiples taxons ids/hits
                        for taxon in temp_subject_taxids:
                            sql = 'INSERT INTO blastnr_taxonomy_%s (' \
                                  'nr_hit_id, hit_number, taxon_id) values (%s, %s, %s)' % (db_name, blast_id, hit_n, taxon)
                            try:
                                server.adaptor.execute(sql)
                            except:
                                print "problem with"
                                print sql

if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb
    import json

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_blast', type=str, help="input interpro xml file", nargs='+')
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-t", '--create_tables', action='store_true', help="Create SQL tables")

    args = parser.parse_args()

    print args.input_blast

    biodb = args.mysql_database

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

    #with open(args.input_taxonomy, 'r') as f:
    #    taxon_id2taxonomy = json.load(f)
    taxon_id2taxonomy = ''

    blastnr2biosql(server,
                    seqfeature_id2locus_tag,
                    locus_tag2genome_taxon_id,
                    protein_id2genome_taxon_id,
                    locus_tag2seqfeature_id,
                    protein_id2seqfeature_id,
                    seqfeature_id2organism,
                    biodb, args.create_tables,
                    *args.input_blast)


