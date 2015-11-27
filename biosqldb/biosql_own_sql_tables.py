#!/usr/bin/python

import manipulate_biosqldb
import gbk2circos

def locus_tag2orthogroup_size(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)
    sql='select locus_tag, count(*) from biosqldb.orthology_detail_%s group by orthogroup'

    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

def superkingdom_table():
    pass



def COG_tables():
    pass

'''


<domain-id>, <genome-name>, <protein-id>,<protein-length>,
<domain-start>, <domain-end>, <COG-id>, <membership-class>,

* Example:

333894695,Alteromonas_SN2_uid67349,333894695,427,1,427,COG0001,0,

CREATE TABLE cog_3014 (domain_id varchar(100),
                        genome_name varchar(200),
                        protein_id INT,
                        protein_length INT,
                        domain_start INT,
                        domain_end INT,
                        COG_id varchar(100),
                        membership_class INT)

# COG	func	name

CREATE table cog_names_2014 (COG_id varchar(100),
                            functon varchar(10),
                             name varchar(200))

create table locus_tag2COG_chlamydia_03_15 (locus_tag varchar(100), COG_id varchar (100))


'''


def create_contig_table(db_name):
    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'CREATE TABLE contigs_%s(accession VARCHAR(100) NOT NULL, ' \
          ' contig_name VARCHAR(100) UNIQUE, ' \
          ' start INT, ' \
          ' end INT)' % db_name


    #server.adaptor.execute(sql,)

    sql2 = 'select accession from bioentry inner join biodatabase on biodatabase.biodatabase_id = bioentry.biodatabase_id' \
           ' where biodatabase.name = "%s"' % db_name
    print sql2
    accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql2,)]



    for accession in accessions:
        record = db.lookup(accession=accession)
        #print record.id
        draft_data = gbk2circos.circos_fasta_draft_misc_features(record)
        if len(draft_data) == 0:
            sql = 'INSERT into contigs_%s (accession, contig_name, start, end) VALUES ("%s", "%s", %s, %s)' % (db_name,
                                                                                                               accession,
                                                                                                                accession,
                                                                                                                0,
                                                                                                                len(record.seq))
            server.adaptor.commit()

            server.adaptor.execute(sql,)
        else:
            for contig in draft_data:
                sql = 'INSERT into contigs_%s (accession, contig_name, start, end) VALUES ("%s", "%s", %s, %s)' % (db_name,
                                                                                                               accession,
                                                                                                               contig[0],
                                                                                                               contig[1],
                                                                                                               contig[2])
                server.adaptor.execute(sql,)
                server.adaptor.commit()


def accession2n_contigs(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select t1.accession, count(*) as n_contigs from bioentry as t1' \
    ' inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id ' \
    ' inner join contigs_chlamydia_03_15 as t3 on t1.accession=t3.accession ' \
    ' where t2.name="%s" group by t3.accession;' % (db_name)

    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


def locus_tag2best_hit(db_name, accession, hit_number=1, rank=False, taxon_name=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    if rank and taxon_name:
        sql_rank_filter = 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.query_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end ' \
                          ' from blastnr.blastnr_hits_%s_%s as t1 ' \
                          ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.%s = "%s";' % (db_name,
                                                                         accession,
                                                                         db_name,
                                                                         accession,
                                                                         db_name,
                                                                         accession,
                                                                         hit_number,
                                                                         rank,
                                                                         taxon_name)
        print sql_rank_filter
        data = server.adaptor.execute_and_fetchall(sql_rank_filter,)

        locus_tag2hit = {}

        for one_hsp in data:
            if not one_hsp[0] in locus_tag2hit:
                locus_tag2hit[one_hsp[0]] = [[ i for i in one_hsp[1:]]]
            else:
                locus_tag2hit[one_hsp[0]].append([ i for i in one_hsp[1:]])

        return locus_tag2hit


    else:
        sql_no_rank_filter = 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.query_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end ' \
                          ' from blastnr.blastnr_hits_%s_%s as t1 ' \
                          ' inner join blastnr.blastnr_hits_taxonomy_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s' % (db_name,
                                                       accession,
                                                       db_name,
                                                       accession,
                                                       db_name,
                                                       accession,
                                                       hit_number)

        print sql_no_rank_filter
        data = server.adaptor.execute_and_fetchall(sql_no_rank_filter,)

        locus_tag2hit = {}

        for one_hsp in data:
            if not one_hsp[0] in locus_tag2hit:
                locus_tag2hit[one_hsp[0]] = [list(one_hsp[1:])]
            else:
                locus_tag2hit[one_hsp[0]].append(list(one_hsp[1:]))

        return locus_tag2hit

def locus_tag2best_hit_n_taxon_ids(db_name, accession):
    server, db = manipulate_biosqldb.load_db(db_name)
    sql = 'select t3.locus_tag, count(*) as n_taxons from blastnr.blastnr_hits_taxonomy_%s_%s as t1 ' \
          ' inner join blastnr.blastnr_hits_%s_%s as t3 on t1.nr_hit_id=t3.nr_hit_id ' \
          ' where t3.hit_number=1 group by t3.nr_hit_id' % (db_name, accession, db_name, accession)
    print sql
    all_blastnr_taxons = server.adaptor.execute_and_fetchall(sql,)



def taxonomical_form(db_name, hit_number=1):

    '''
    contruct of an html form with accessions, superkingdom and phylum selection
    possibility to nest more levels

    '''

    server, db = manipulate_biosqldb.load_db(db_name)

    sql1 = 'select accession, bioentry.description from biosqldb.bioentry' \
           ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id' \
           ' where biodatabase.name = "%s"' % db_name

    accessions = server.adaptor.execute_and_fetchall(sql1,)

    #accessions = [accession + ('-',) for accession in accessions]


    form = '<p><label for="id_genome_taxonomy">Genome:</label><select id="genome" name="Genome">\n'
    form+='    <option value="">--</option>\n'
    for accession in accessions:
        form+='<option value="%s">%s</option>\n' % (accession[0], accession[1])

    form+='</select></p>\n'
    #print form

    # pour chaque accession, faire un choix avec les superkingdom

    form+='<p><label for="id_superkingdom_taxonomy">Superkingdom:</label><select id="superkingdom" name="Superkingdom">\n'
    form+='    <option value="">--</option>\n'
    for accession in accessions:

        sql2 = 'select t3.superkingdom ' \
              ' from blastnr.blastnr_hits_%s_%s as t1 ' \
              ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
              ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
              '  where t1.hit_number=%s group by t3.superkingdom' % (db_name, accession[0],db_name, accession[0], hit_number)

        #print sql2

        superkinkdoms = [i[0] for i in server.adaptor.execute_and_fetchall(sql2,)]

        for superkinkdom in superkinkdoms:
            if superkinkdom == "-":
                #print superkinkdoms, accession
                form+='    <option value="%s_%s" class="%s">%s</option>\n' % (accession[0], 'unclassified', accession[0], 'unclassified')
            else:
                form+='    <option value="%s_%s" class="%s">%s</option>\n' % (accession[0], superkinkdom, accession[0], superkinkdom)
    form+='</select></p>\n'



    form+='<p><label for="id_phylum_taxonomy">Phylum:</label><select id="phylum" name="Phylum">\n'
    form+='    <option value="">--</option>\n'
    for accession in accessions:

        sql3 = 'select distinct t3.superkingdom, t3.phylum ' \
              ' from blastnr.blastnr_hits_%s_%s as t1 ' \
              ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
              ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
              ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
              ' where t1.hit_number=%s' % (db_name, accession[0],db_name, accession[0], db_name, accession[0], hit_number)
        #print sql3
        phylums = server.adaptor.execute_and_fetchall(sql3,)

        #print sql3

        for phylum in phylums:
            #print accession, phylum
            if phylum[0] == '-' and phylum[1] != '-':
                #print 'cas 1', accession
                form+='    <option value="%s_%s_%s" class="%s_%s">%s</option>\n' % (accession[0], 'unclassified', phylum[1], accession[0],'unclassified', phylum[1])
            elif phylum[0] == '-' and phylum[1] == '-':
                #print 'cas 2', accession

                form+='    <option value="%s_%s_%s" class="%s_%s">%s</option>\n' % (accession[0], 'unclassified', 'unclassified', accession[0],'unclassified', 'unclassified')
            elif phylum[0] != '-' and phylum[1] == '-':
                #print 'cas 3', accession
                form+='    <option value="%s_%s_%s" class="%s_%s">%s</option>\n' % (accession[0], phylum[0], 'unclassified', accession[0], phylum[0], 'unclassified')



            else:
                form+='    <option value="%s_%s_%s" class="%s_%s">%s</option>\n' % (accession[0], phylum[0], phylum[1], accession[0], phylum[0], phylum[1])
    form+='</select></p>\n'
    return form

def clean_multispecies_blastnr_record(db_name, create_new_sql_tables = False):


    server, db = manipulate_biosqldb.load_db(db_name)
    sql = 'select accession from bioentry' \
      ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id' \
      ' and biodatabase.name = "%s"' % db_name

    all_accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    #all_accessions = ['Rhab']
    if create_new_sql_tables:
        for accession in all_accessions:
            sql_blast_taxonomy = 'CREATE TABLE blastnr.blastnr_hits_taxonomy_filtered_%s_%s (nr_hit_id INT, ' \
                         'subject_taxon_id int)' % (db_name, accession)

            sql_blast_taxonomy2 = 'ALTER TABLE blastnr.blastnr_hits_taxonomy_filtered_%s_%s ADD CONSTRAINT fk_blast_hit_id_filter_%s ' \
                                  'FOREIGN KEY (nr_hit_id) REFERENCES blastnr.blastnr_hits_%s_%s(nr_hit_id);' % (db_name, accession, accession, db_name, accession)
            print sql_blast_taxonomy
            print sql_blast_taxonomy2
            server.adaptor.execute(sql_blast_taxonomy,)
            server.adaptor.execute(sql_blast_taxonomy2,)
            server.adaptor.commit()


    total_ambiguous_taxons = 0
    for accession in all_accessions:
        #print accession, total_ambiguous_taxons
        sql = 'select nr_hit_id from blastnr.blastnr_hits_%s_%s as t1' % (db_name, accession) #  where hit_number=1

        all_blastnr_hit_ids = server.adaptor.execute_and_fetchall(sql,)

        more_than_one_genus_path = 0
        for blast_id in all_blastnr_hit_ids:
            blast_id = int(blast_id[0])

            sql = 'select  t3.locus_tag, t1.nr_hit_id, t1.subject_taxon_id,  t2.phylum, t2.order, t2.family, t2.genus, t2.species, t2.superkingdom, t3.hit_number, t3.subject_title from blastnr.blastnr_hits_taxonomy_%s_%s as t1 ' \
                  ' inner join blastnr.blastnr_hits_%s_%s as t3 on t1.nr_hit_id=t3.nr_hit_id ' \
                  ' inner join blastnr.blastnr_taxonomy as t2 on t1.subject_taxon_id=t2.taxon_id ' \
                  ' where t1.nr_hit_id=%s' % (db_name, accession, db_name, accession, blast_id)
            one_blastnr_taxons = server.adaptor.execute_and_fetchall(sql,)
            #print 'blast_id', blast_id, 'n taxons', len(one_blastnr_taxons)

            # build dictionnary
            one_nr_hit_taxons = {}
            for i in one_blastnr_taxons:
                one_nr_hit_taxons[i[2]] = {
                    'query_locus_tag': i[0],
                    'nr_hit_id': i[1],
                    'phylum': i[3],
                    'order': i[4],
                    'family': i[5],
                    'genus': i[6],
                    'species': i[7],
                    'superkingdom': i[8],
                    'hit_number': i[9],
                    'description': i[10]
                }

            '''
            best_classified_path_length = 0
            for one_taxon in one_blastnr_taxons:

                # flag to check if redundancy was not already observed in the taxon group
                temp_count = 0
                for one_rank in one_taxon[3:]:
                    if one_rank != '-':
                        temp_count += 1
                if temp_count > best_classified_path_length:
                    best_classified_path_length = temp_count

            # equally lon path extraction
            '''
            unique_taxon_paths = {}
            poor_path = {}
            # iter dictionnary
            for one_subject_taxon_id in one_nr_hit_taxons:

                # removing data from the phage (small contig identical to the phage)
                phage = ['Rhab_01766','Rhab_01767','Rhab_01768','Rhab_01769','Rhab_01770','Rhab_01771','Rhab_01772']
                #print one_nr_hit_taxons[one_subject_taxon_id]
                if one_nr_hit_taxons[one_subject_taxon_id]['query_locus_tag'] in phage:
                    break

                # check if phylum, family, genus, species are not null
                temp_count = 0
                for one_rank in ['phylum', 'order', 'family', 'genus', 'species', 'superkingdom']:
                    if one_nr_hit_taxons[one_subject_taxon_id][one_rank] != '-':
                        temp_count += 1
                one_nr_hit_taxons[one_subject_taxon_id]['path_completeness'] = temp_count


                # check if an idenical subject taxon path is already present in equall_taxon_path
                # different taxon id can have identical phylum-genus path: i.e strains of the same species, subspecies,..
                if one_subject_taxon_id not in unique_taxon_paths:

                    counted = False
                    # check if taxon path with same family and same genus is present
                    # check if taxon path with same family and no genus ('-') is present
                    ### => replace it if genus exist in the new path
                    for unique_taxon_path in unique_taxon_paths.keys():
                        if (unique_taxon_paths[unique_taxon_path]['family'] == one_nr_hit_taxons[one_subject_taxon_id]['family']
                            and unique_taxon_paths[unique_taxon_path]['genus'] == one_nr_hit_taxons[one_subject_taxon_id]['genus']) \
                                or (unique_taxon_paths[unique_taxon_path]['family'] == one_nr_hit_taxons[one_subject_taxon_id]['family']
                                    and one_nr_hit_taxons[one_subject_taxon_id]['genus'] == '-'):
                            unique_taxon_paths[unique_taxon_path]['taxon_path_frequency'] += 1
                            counted = True
                            # one match, break the loop
                            break
                        # if same family as the reference but the reference has no genus and the new taxon has one
                        elif (unique_taxon_paths[unique_taxon_path]['family'] == one_nr_hit_taxons[one_subject_taxon_id]['family']
                                and unique_taxon_paths[unique_taxon_path]['genus'] == '-'):

                            # retrieve the count of the uncomplete path, replace it and delete it
                            unique_taxon_paths[one_subject_taxon_id] = one_nr_hit_taxons[one_subject_taxon_id]
                            unique_taxon_paths[one_subject_taxon_id]['taxon_path_frequency'] = unique_taxon_paths[unique_taxon_path]['taxon_path_frequency'] + 1
                            del unique_taxon_paths[unique_taxon_path]
                            counted = True
                            # one match break the loop
                            break
                    # if no match with any unique path, create a new unique entry with count 1 (new family / genus combination)
                    if not counted:
                        unique_taxon_paths[one_subject_taxon_id] = one_nr_hit_taxons[one_subject_taxon_id]
                        unique_taxon_paths[one_subject_taxon_id]['taxon_path_frequency'] = 1
                # if two time identical taxon (can not happen, oder?)
                else:
                    unique_taxon_paths[one_subject_taxon_id]['taxon_path_frequency']  += 1

            '''
            # count the most frequent paths and keep all equally frequent path
            if len(equall_taxon_path) > 1:
                max_frequency = 0
                for i in equall_taxon_path:
                    if equall_taxon_path[i][0] > max_frequency:
                        max_frequency = equall_taxon_path[i][0]

                keep_frequency = []
                for i in equall_taxon_path:
                    if equall_taxon_path[i][0] == max_frequency:
                        keep_frequency.append(equall_taxon_path[i])


                for i in keep_frequency:
                    print 'frequency:', i[0], 'path length', i[2], i[1][2:]
            '''



            poor_path = {}
            ok_path = {}
            # if with have path with different family/genus
            # check if all some are not complete (if only 2 rank/5 are recorded) ('chlamydiae', 'chlamydiales', '-', '-', '-', '-')
            if len(unique_taxon_paths) > 1:

                for one_path_taxid in unique_taxon_paths:
                    if unique_taxon_paths[one_path_taxid]['path_completeness'] < 4:
                        #print 'poor path!', unique_taxon_paths[one_path_taxid]
                        poor_path[one_path_taxid] = unique_taxon_paths[one_path_taxid]
                    else:
                        ok_path[one_path_taxid] = unique_taxon_paths[one_path_taxid]
                # if we have at least one good path

                if len(ok_path) > 1:

                    # hit associated with different families/genus
                    more_than_one_genus_path += 1
                    for taxon_id in ok_path:

                        sql = 'INSERT INTO blastnr.blastnr_hits_taxonomy_filtered_%s_%s (' \
                        'nr_hit_id, subject_taxon_id) values (%s, %s)'  % (db_name,
                                                                           accession,
                                                                           ok_path[taxon_id]['nr_hit_id'],
                                                                           taxon_id)

                        server.adaptor.execute(sql,)
                        server.adaptor.commit()
                        '''
                        print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (len(ok_path),
                                    ok_path[taxon_id]['nr_hit_id'],
                                    ok_path[taxon_id]['hit_number'],
                                    ok_path[taxon_id]['description'],
                                    ok_path[taxon_id]['query_locus_tag'],
                                    taxon_id,
                                    ok_path[taxon_id]['superkingdom'],
                                    ok_path[taxon_id]['phylum'],
                                    ok_path[taxon_id]['order'],
                                    ok_path[taxon_id]['family'],
                                    ok_path[taxon_id]['genus'])
                        '''

                # if we do not have any good path
                elif len(ok_path) == 0 and len(poor_path) == 0:
                    #print 'No path: Impossible!'
                    import sys
                    sys.exit()
                elif len(ok_path) == 1:

                    for taxon_id in ok_path:

                        sql = 'INSERT INTO blastnr.blastnr_hits_taxonomy_filtered_%s_%s (' \
                        'nr_hit_id, subject_taxon_id) values (%s, %s)'  % (db_name,
                                                                           accession,
                                                                           ok_path[taxon_id]['nr_hit_id'],
                                                                           taxon_id)

                        server.adaptor.execute(sql,)
                        server.adaptor.commit()
                else:
                    # no good path
                    for taxon_id in poor_path:

                        sql = 'INSERT INTO blastnr.blastnr_hits_taxonomy_filtered_%s_%s (' \
                        'nr_hit_id, subject_taxon_id) values (%s, %s)'  % (db_name,
                                                                           accession,
                                                                           poor_path[taxon_id]['nr_hit_id'],
                                                                           taxon_id)

                        server.adaptor.execute(sql,)
                        server.adaptor.commit()

            # only one path
            else:
                for taxon_id in unique_taxon_paths:

                    sql = 'INSERT INTO blastnr.blastnr_hits_taxonomy_filtered_%s_%s (' \
                        'nr_hit_id, subject_taxon_id) values (%s, %s)'  % (db_name,
                                                                           accession,
                                                                           unique_taxon_paths[taxon_id]['nr_hit_id'],
                                                                           taxon_id)

                    server.adaptor.execute(sql,)
                    server.adaptor.commit()


    #print "total_ambiguous_taxons", total_ambiguous_taxons

    '''
    nr_hit_id |
    subject_taxon_id |
    taxon_id |

    no_rank            |
    superkingdom |
    kingdom |
    subkingdom |
    superphylum |
    phylum         |
    subphylum |
    superclass |
    class               |
    subclass |
    superorder |
    order           |
    suborder |
    superfamily |
    family           |
    subfamily |
    genus         |
    subgenus |
    species |
    species_subgroup |
    species_group |
    subspecies |
    tribe |
    infraorder |
    subtribe |
    forma |
    infraclass |
    varietas |
    parvorder |
    '''



def get_best_hit_excluding_one_family():
    sql='select t3.subject_taxon_id,t4.no_rank, t4.phylum, t4.order, t4.family, B.* from (select * from biosqldb.orthology_detail_chlamydia_03_15 as t1 where t1.accession="NC_015713") A left join (select t2.nr_hit_id,t2.locus_tag,t2.subject_kingdom,t2.subject_accession from blastnr.blastnr_hits_chlamydia_03_15_NC_015713 as t2) B on A.locus_tag=B.locus_tag inner join blastnr.blastnr_hits_taxonomy_chlamydia_03_15_NC_015713 as t3 on B.nr_hit_id=t3.nr_hit_id inner join blastnr.blastnr_taxonomy as t4 on t3.subject_taxon_id=t4.taxon_id where t4.family !="Simkaniaceae";'

def locus_tag2n_nr_hits(db_name, genome_accession, exclude_family = False):
    print 'exclude family', exclude_family
    server, db = manipulate_biosqldb.load_db(db_name)
    if not exclude_family:
        sql = 'select  A.locus_tag, B.n_hits from (select * from biosqldb.orthology_detail_%s as t1 where t1.accession="%s") A ' \
              ' left join (select t2.nr_hit_id,t2.locus_tag,t2.subject_kingdom,t2.subject_accession,count(t2.locus_tag) as n_hits ' \
              ' from blastnr.blastnr_hits_%s_%s as t2 ' \
              ' group by locus_tag) B on A.locus_tag=B.locus_tag;' %(db_name, genome_accession, db_name, genome_accession)
    else:
        sql = 'select  A.locus_tag, B.n_hits from (select * from biosqldb.orthology_detail_%s as t1 where t1.accession="%s") A ' \
              ' left join (select t2.nr_hit_id,t2.locus_tag,t2.subject_kingdom,t2.subject_accession, count(t2.locus_tag) as n_hits ' \
              ' from blastnr.blastnr_hits_%s_%s as t2 inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t3 on t3.nr_hit_id=t2.nr_hit_id' \
              ' left join blastnr.blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id where t4.family!="%s"' \
              ' group by locus_tag) B on A.locus_tag=B.locus_tag;' % (db_name,
                                                                     genome_accession,
                                                                     db_name,
                                                                     genome_accession,
                                                                     db_name,
                                                                     genome_accession,
                                                                     exclude_family)


    print sql
    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

def best_hit_classification(db_name, accession):
    sql = 'select A.*, t4.superkingdom, t4.kingdom, t4.phylum, t4.order, t4.family, t4.genus, t4.species' \
          ' from (select locus_tag, TM, SP, gene, product from biosqldb.orthology_detail_%s as t1 where t1.accession="%s" group by locus_tag) A' \
          ' left join (select t2.nr_hit_id, locus_tag from blastnr.blastnr_hits_%s_%s as t2  where hit_number=1) B on A.locus_tag=B.locus_tag' \
          ' left join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
          ' left join blastnr.blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id;' % (db_name,
                                                                                                                        accession,
                                                                                                                        db_name,
                                                                                                                        accession,                                                                                                              db_name,
                                                                                                                        accession)

    import pandas
    server, db = manipulate_biosqldb.load_db(db_name)

    data = [i for i in server.adaptor.execute_and_fetchall(sql)]
    table = pandas.DataFrame(data, columns=['locus_tag', 'TM', 'SP', 'gene', 'product', 'superkingdom','kingdom','phylum','order','family','genus','species'])
    return table

def best_hit_phylum_and_protein_length(db_name, accession):
    sql = 'select A.locus_tag, char_length(A.translation), t4.superkingdom' \
          ' from (select locus_tag, translation from biosqldb.orthology_detail_%s as t1 where t1.accession="%s" group by locus_tag) A' \
          ' left join (select t2.nr_hit_id, locus_tag from blastnr.blastnr_hits_%s_%s as t2  where hit_number=1) B on A.locus_tag=B.locus_tag' \
          ' left join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
          ' left join blastnr.blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id;' % (db_name,
                                                                                                                        accession,
                                                                                                                        db_name,
                                                                                                                        accession,                                                                                                              db_name,
                                                                                                                        accession)

    import pandas
    server, db = manipulate_biosqldb.load_db(db_name)

    data = [i for i in server.adaptor.execute_and_fetchall(sql)]
    table = pandas.DataFrame(data, columns=['locus_tag', 'protein_length', 'superkingdom'])
    return table



def locus_tag2n_blast_superkingdom(db_name, accession, superkingdom="Bacteria", exclude_family=False):
    if not exclude_family:
        sql = 'select A.locus_tag, count(*)' \
              ' from (select locus_tag, translation from biosqldb.orthology_detail_%s as t1 where t1.accession="%s" group by locus_tag) A' \
              ' left join (select t2.nr_hit_id, locus_tag from blastnr.blastnr_hits_%s_%s as t2) B on A.locus_tag=B.locus_tag' \
              ' left join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
              ' left join blastnr.blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
              ' where t4.superkingdom="%s" group by locus_tag;' % (db_name,
                                                                        accession,
                                                                        db_name,
                                                                        accession,                                                                                                              db_name,
                                                                        accession,
                                                                        superkingdom)
    else:
        sql = 'select A.locus_tag, count(*)' \
              ' from (select locus_tag, translation from biosqldb.orthology_detail_%s as t1 where t1.accession="%s" group by locus_tag) A' \
              ' left join (select t2.nr_hit_id, locus_tag from blastnr.blastnr_hits_%s_%s as t2) B on A.locus_tag=B.locus_tag' \
              ' left join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
              ' left join blastnr.blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
              ' where t4.superkingdom="%s" and t4.family!="%s" group by locus_tag;' % (db_name,
                                                                        accession,
                                                                        db_name,
                                                                        accession,                                                                                                              db_name,
                                                                        accession,
                                                                        superkingdom,
                                                                        exclude_family)



    server, db = manipulate_biosqldb.load_db(db_name)
    data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    return data

def locus_tag2n_blast_bacterial_phylum(db_name, accession, phylum="Chlamydiae", reverse=False, exclude_family = False):
    if not reverse:
        if not exclude_family:
            sql = 'select A.locus_tag, count(*)' \
                  ' from (select locus_tag, translation from biosqldb.orthology_detail_%s as t1 where t1.accession="%s" group by locus_tag) A' \
                  ' left join (select t2.nr_hit_id, locus_tag from blastnr.blastnr_hits_%s_%s as t2) B on A.locus_tag=B.locus_tag' \
                  ' left join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
                  ' left join blastnr.blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
                  ' where t4.superkingdom="Bacteria" and t4.phylum="%s" group by locus_tag;' % (db_name,
                                                                            accession,
                                                                            db_name,
                                                                            accession,                                                                                                              db_name,
                                                                            accession,
                                                                            phylum)
        else:
            sql = 'select A.locus_tag, count(*)' \
                  ' from (select locus_tag, translation from biosqldb.orthology_detail_%s as t1 where t1.accession="%s" group by locus_tag) A' \
                  ' left join (select t2.nr_hit_id, locus_tag from blastnr.blastnr_hits_%s_%s as t2) B on A.locus_tag=B.locus_tag' \
                  ' left join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
                  ' left join blastnr.blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
                  ' where t4.superkingdom="Bacteria" and t4.phylum="%s" and t4.family!="%s" group by locus_tag;' % (db_name,
                                                                            accession,
                                                                            db_name,
                                                                            accession,                                                                                                              db_name,
                                                                            accession,
                                                                            phylum,
                                                                            exclude_family)
    else:
        if not exclude_family:
            sql = 'select A.locus_tag, count(*)' \
                  ' from (select locus_tag, translation from biosqldb.orthology_detail_%s as t1 where t1.accession="%s" group by locus_tag) A' \
                  ' left join (select t2.nr_hit_id, locus_tag from blastnr.blastnr_hits_%s_%s as t2) B on A.locus_tag=B.locus_tag' \
                  ' left join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
                  ' left join blastnr.blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
                  ' where t4.superkingdom="Bacteria" and t4.phylum !="%s" group by locus_tag;' % (db_name,
                                                                            accession,
                                                                            db_name,
                                                                            accession,                                                                                                              db_name,
                                                                            accession,
                                                                            phylum)
        else:
            sql = 'select A.locus_tag, count(*)' \
                  ' from (select locus_tag, translation from biosqldb.orthology_detail_%s as t1 where t1.accession="%s" group by locus_tag) A' \
                  ' left join (select t2.nr_hit_id, locus_tag from blastnr.blastnr_hits_%s_%s as t2) B on A.locus_tag=B.locus_tag' \
                  ' left join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
                  ' left join blastnr.blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
                  ' where t4.superkingdom="Bacteria" and t4.phylum !="%s" and t4.family!="%s" group by locus_tag;' % (db_name,
                                                                            accession,
                                                                            db_name,
                                                                            accession,                                                                                                              db_name,
                                                                            accession,
                                                                            phylum, exclude_family)

    server, db = manipulate_biosqldb.load_db(db_name)
    data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    return data







def locus_tag2orthogroup(db_name):
    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select locus_tag, orthogroup from orthology_detail_%s group by locus_tag' % db_name

    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


def locus_tag2presence_in_n_genomes(db_name):
    '''
    return a dictionnary of all locus tag and the number of genomes in which they have one or multiple homolog(s)
    '''

    import mysqldb_load_mcl_output
    server, db = manipulate_biosqldb.load_db(db_name)
    orthogroup2family_size = mysqldb_load_mcl_output.get_family_size(server, db_name)
    locus_tag2orthogroup_dico = locus_tag2orthogroup(db_name)

    locus_tag2presence_absence = {}
    for locus in locus_tag2orthogroup_dico:
        locus_tag2presence_absence[locus] = orthogroup2family_size[locus_tag2orthogroup_dico[locus]]
    return locus_tag2presence_absence

def locus_tag2n_paralogs(db_name, genome_accession):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select orthogroup, count(*) as paralogs ' \
          ' from biosqldb.orthology_detail_%s ' \
          ' where accession="%s" group by orthogroup order by paralogs DESC;' % (db_name, genome_accession)

    return server.adaptor.execute_and_fetchall(sql,)

def orthogroup2gene(db_name, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    if not accession:
        sql =  'select orthogroup, gene from orthology_detail_%s' % db_name
    else:
        sql = 'select orthogroup, gene from orthology_detail_%s where accession="%s"' % (db_name, accession)


    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2gene = {}

    for row in data:
        if row[0] not in ortho2gene:
            ortho2gene[row[0]] = {}
            ortho2gene[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2gene[row[0]]:
                ortho2gene[row[0]][row[1]] += 1
            else:
                ortho2gene[row[0]][row[1]] = 1

    return ortho2gene






def orthogroup2cog(db_name, group_list, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    group_list_form = '"' + '","'.join(group_list) + '"'

    if not accession:
        sql =  'select t1.orthogroup, t2.cog_id from (select * from biosqldb.orthology_detail_%s ' \
               ' where orthogroup in (%s)) t1 left join COG.locus_tag2gi_hit_%s ' \
               'as t2 on t1.locus_tag=t2.locus_tag;' % (db_name, group_list_form, db_name)
    else:
        sql = 'select t1.orthogroup, t2.cog_id from biosqldb.orthology_detail_%s as t1 left join COG.locus_tag2gi_hit_%s ' \
              'as t2 on t1.locus_tag=t2.locus_tag' % (db_name, db_name)

    print 'df', sql
    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2cog = {}

    for row in data:
        if row[0] not in ortho2cog:
            ortho2cog[row[0]] = {}
            ortho2cog[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2cog[row[0]]:
                ortho2cog[row[0]][row[1]] += 1
            else:
                ortho2cog[row[0]][row[1]] = 1

    return ortho2cog

def orthogroup2pfam(db_name, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    if not accession:
        sql =  'select orthogroup, signature_accession from biosqldb.interpro_%s ' \
               ' where analysis="Pfam"' % (db_name, )
    else:
        sql = 'select t1.orthogroup, t2.cog_id from biosqldb.orthology_detail_%s as t1 left join (COG.locus_tag2gi_hit_%s ' \
              'as t2 on t1.locus_tag=t2.locus_tag' % (db_name, db_name)

    print 'df', sql
    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2pfam = {}

    for row in data:
        if row[0] not in ortho2pfam:
            ortho2pfam[row[0]] = {}
            ortho2pfam[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2pfam[row[0]]:
                ortho2pfam[row[0]][row[1]] += 1
            else:
                ortho2pfam[row[0]][row[1]] = 1

    return ortho2pfam


def orthogroup2interpro(db_name, group_list, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    group_list_form = '"' + '","'.join(group_list) + '"'

    if not accession:
        sql =  'select t1.orthogroup, t2.cog_id from (select * from biosqldb.orthology_detail_%s ' \
               ' where orthogroup in (%s)) t1 left join COG.locus_tag2gi_hit_%s ' \
               'as t2 on t1.locus_tag=t2.locus_tag;' % (db_name, group_list_form, db_name)
    else:
        sql = 'select t1.orthogroup, t2.cog_id from biosqldb.orthology_detail_%s as t1 left join COG.locus_tag2gi_hit_%s ' \
              'as t2 on t1.locus_tag=t2.locus_tag' % (db_name, db_name)

    print 'df', sql
    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2cog = {}

    for row in data:
        if row[0] not in ortho2cog:
            ortho2cog[row[0]] = {}
            ortho2cog[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2cog[row[0]]:
                ortho2cog[row[0]][row[1]] += 1
            else:
                ortho2cog[row[0]][row[1]] = 1

    return ortho2cog


def pfam2description(db_name, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)


    if not accession:
        sql =  'select signature_accession, signature_description from biosqldb.interpro_%s ' \
               ' where analysis="Pfam" group by signature_accession' % (db_name)
    else:
        sql = 'select t1.orthogroup, t2.cog_id from biosqldb.orthology_detail_%s as t1 left join COG.locus_tag2gi_hit_%s ' \
              'as t2 on t1.locus_tag=t2.locus_tag' % (db_name, db_name)


    data = server.adaptor.execute_and_fetchall(sql,)



    return manipulate_biosqldb.to_dict(data)









def orthogroup2product(db_name, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    if not accession:
        sql =  'select orthogroup, product from orthology_detail_%s' % db_name
    else:
        sql = 'select orthogroup, product from orthology_detail_%s where accession="%s"' % (db_name, accession)


    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2product = {}

    for row in data:
        if row[0] not in ortho2product:
            ortho2product[row[0]] = {}
            ortho2product[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2product[row[0]]:
                ortho2product[row[0]][row[1]] += 1
            else:
                ortho2product[row[0]][row[1]] = 1

    return ortho2product


'''
print 'tata'
locus_tag2best_hit_dico = locus_tag2best_hit("chlamydia_03_15", "Rhab", hit_number=1, rank=False, taxon_name=False)
print
print locus_tag2best_hit_dico.keys()[0]
print locus_tag2best_hit_dico[locus_tag2best_hit_dico.keys()[0]]

locus_tag2best_hit_euk = locus_tag2best_hit("chlamydia_03_15", "Rhab", hit_number=1, rank="superkingdom", taxon_name="Eukaryota")
print locus_tag2best_hit_euk
'''

def calculate_average_protein_identity(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select taxon_id from bioentry' \
          ' inner join biodatabase on bioentry.biodatabase_id=biodatabase.biodatabase_id where biodatabase.name="%s" group by taxon_id' % db_name

    all_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


    average_id = {}
    for i_taxon in range(0,len(all_taxons)):
        print 'i taxon', i_taxon
        for y_taxon in range(i_taxon+1, len(all_taxons)):
            shared_groups_sql = 'select orthogroup from orthology_%s where `%s` =1 and `%s`=1' % (db_name, all_taxons[i_taxon], all_taxons[y_taxon])
            all_groups = [i[0] for i in server.adaptor.execute_and_fetchall(shared_groups_sql,)]
            identity_values = []
            print 'y taxon', y_taxon
            for group in all_groups:
                sql_ref = 'select locus_tag from orthology_detail_%s where taxon_id=%s and orthogroup="%s"' % (db_name, all_taxons[i_taxon], group)
                sql_query = 'select locus_tag from orthology_detail_%s where taxon_id=%s and orthogroup="%s"' % (db_name, all_taxons[y_taxon], group)
                locus_ref = server.adaptor.execute_and_fetchall(sql_ref,)[0][0]
                locus_query = server.adaptor.execute_and_fetchall(sql_query,)[0][0]

                sql_id = 'select `%s` from orth_%s.%s where locus_tag="%s"' % (locus_ref, db_name, group, locus_query)

                id_value = server.adaptor.execute_and_fetchall(sql_id,)[0][0]

                identity_values.append(id_value)
            one_average_id = sum(identity_values) / float(len(identity_values))
            if not all_taxons[i_taxon] in average_id:
                average_id[all_taxons[i_taxon]] = {}
            average_id[all_taxons[i_taxon]][all_taxons[y_taxon]] = [one_average_id, len(identity_values)]
            print "one_average_id",all_taxons[i_taxon],all_taxons[y_taxon], one_average_id, len(identity_values)
    import json
    with open('genome_identity_dico.json', 'wb') as fp:
        json.dump(average_id, fp)


def get_cooccurring_groups(db_name, accession1, accession2, windows_size=10, min_number_cooccurring=4, exclude=[]):
    pass
    server, db = manipulate_biosqldb.load_db(db_name)

    reference = db.lookup(accession=accession1)
    query = db.lookup(accession=accession2)

    def get_orthogroups(record):
        groups = []
        for i in record.features:
            if i.type == 'CDS':
                try:
                    groups.append(i.qualifiers['orthogroup'][0])
                except:
                    pass

        return groups

    reference_groups = get_orthogroups(reference)
    query_groups = get_orthogroups(query)


    conserved_regions = []
    for start in range(0, len(reference_groups)):
        ref_window = reference_groups[start:start+windows_size]
        for start_query in range(0, len(query_groups)):
            query_window = query_groups[start_query:start_query+windows_size]
            #print query_window, ref_window
            identical_groups = list(set(ref_window).intersection(query_window))
            if len(identical_groups) >= min_number_cooccurring:
                for i in identical_groups:
                    if i in exclude:
                        break
                else:
                    if identical_groups not in conserved_regions:
                        conserved_regions.append(identical_groups)

    return conserved_regions

def orthogroup2protein_id_list(db_name):
    server, db = manipulate_biosqldb.load_db(db_name)

    sql1 = 'select orthogroup, protein_id, locus_tag from orthology_detail_%s;' % db_name

    data1 = server.adaptor.execute_and_fetchall(sql1)

    orthogroup2protein_id_list_dico = {}
    for i in data1:
        orthogroup = i[0]
        if i[1] == '-':
            protein_id = i[2]
        else:
            protein_id = i[1]
        if not orthogroup in orthogroup2protein_id_list_dico:
            orthogroup2protein_id_list_dico[orthogroup] = [protein_id]
        else:
            orthogroup2protein_id_list_dico[orthogroup].append(protein_id)

    return orthogroup2protein_id_list_dico





if __name__ == '__main__':
    #print taxonomical_form('chlamydia_03_15')
    #locus_tag2best_hit_n_taxon_ids('chlamydia_03_15', 'Rhab')
    #clean_multispecies_blastnr_record('chlamydia_03_15', create_new_sql_tables=False)

    #locus_tag2n_nr_hits('chlamydia_03_15', 'Rhab')
    #locus_tag2n_blast_bacteria('chlamydia_03_15', 'Rhab')
    #calculate_average_protein_identity('chlamydia_03_15')
    #
    '''
    ribosomal_proteins = get_cooccurring_groups('chlamydia_03_15', 'Rhab', 'CP001928', 10, 9)
    merged_ribosomal_proteins = []
    for i in ribosomal_proteins:
        merged_ribosomal_proteins += i

    test = ribosomal_proteins = get_cooccurring_groups('chlamydia_03_15', 'Rhab', 'CP001928', 10, 7, merged_ribosomal_proteins)
    print len(test)
    for i in test:
        print i
    '''
    wcw_syntheny = get_cooccurring_groups('chlamydia_03_15', 'Rhab', 'CP001928',3,3)
    simkania_syntheny = get_cooccurring_groups('chlamydia_03_15', 'Rhab', 'NC_015713',3,3)
    trachomatis_syntheny = get_cooccurring_groups('chlamydia_03_15', 'Rhab', 'AE001273',3,3)
    para_proto_syntheny = get_cooccurring_groups('chlamydia_03_15', 'NC_005861', 'KNic',3,3)
    tracho_muridarum = get_cooccurring_groups('chlamydia_03_15', 'NC_002620', 'AE001273',3,3)
    print len(wcw_syntheny), len(simkania_syntheny), len(trachomatis_syntheny)
    print len(para_proto_syntheny), len(tracho_muridarum)


def _chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]


def update_interpro_table(biodb, locus_list, locus2ortho, i):
    server, db = manipulate_biosqldb.load_db(biodb)
    y = 0
    for locus in locus_list:
        y+=1
        print i,y, len(locus_list), locus
        sql = 'UPDATE interpro_%s SET orthogroup="%s" WHERE  locus_tag="%s";' % (biodb, locus2ortho[locus], locus)
        server.adaptor.execute(sql,)
        server.commit()

def add_orthogroup_to_interpro_table(biodb_name):

    '''
    Drop table interpro, puis on la reconstruit
    '''

    locus2ortho = locus_tag2orthogroup(biodb_name)

    import numpy
    from multiprocessing import Process
    import time

    server, db = manipulate_biosqldb.load_db(biodb_name)
    #sql = 'ALTER TABLE interpro_%s ADD orthogroup VARCHAR(100);' % biodb_name
    sql = 'select * from interpro_%s;' % biodb_name

    interpro_data = server.adaptor.execute_and_fetchall(sql,)

    try:
        sql = 'drop table interpro_%s' % biodb_name
        server.adaptor.execute(sql,)
        server.adaptor.commit()
    except:
        pass

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
          ' pathways varchar(10000),' \
          ' orthogroup varchar(100))' % biodb_name
    server.adaptor.execute(sql)
    i=0
    for one_row in interpro_data:
        print i
        i+=1
        sql = 'INSERT INTO interpro_%s(accession, locus_tag, organism, taxon_id,' \
              ' sequence_length, analysis, signature_accession, signature_description, start, ' \
              ' stop, score, interpro_accession, interpro_description, GO_terms, pathways, orthogroup) ' \
              ' values ("%s", "%s", "%s", %s, %s, "%s", "%s", "%s", %s, %s, "%s", "%s", "%s", "%s", "%s", "%s");' % (biodb_name,
                                                                                             one_row[0],
                                                                                             one_row[1],
                                                                                             one_row[2],
                                                                                             one_row[3],
                                                                                       one_row[4],
                                                                                       one_row[5],
                                                                                       one_row[6],
                                                                                       one_row[7],
                                                                                       one_row[8],
                                                                                       one_row[9],
                                                                                       one_row[10],
                                                                                       one_row[11],
                                                                                       one_row[12],
                                                                                       one_row[13],
                                                                                       one_row[14],
                                                                                       locus2ortho[one_row[1]])
        #try:

        server.adaptor.execute(sql)
        server.adaptor.commit()
        #except:
        #    pass

    '''
    locus2ortho = locus_tag2orthogroup(biodb_name)
    n_cpu = 8
    n_poc_per_list = int(numpy.ceil(len(locus2ortho)/float(n_cpu)))
    query_lists = _chunks(locus2ortho.keys(), n_poc_per_list)
    print len(query_lists)
    time.sleep(5)
    procs = []
    i = 0
    for one_list in query_lists:
        i+= 1
        proc = Process(target=update_interpro_table, args=(biodb_name, one_list, locus2ortho, i))#, out_q))
        procs.append(proc)
        proc.start()

    time.sleep(5)
    for proc in procs:
        proc.join()
    '''



    '''

    import json
    import re
    server, db = manipulate_biosqldb.load_db('chlamydia_03_15')
    taxon_id2genome_description = manipulate_biosqldb.taxon_id2genome_description(server, 'chlamydia_03_15')

    for i, accession in enumerate(taxon_id2genome_description):
        #print i, accession
        description = taxon_id2genome_description[accession]
        description = re.sub(", complete genome\.", "", description)
        description = re.sub(", complete genome", "", description)
        description = re.sub(", complete sequence\.", "", description)
        description = re.sub("strain ", "", description)
        description = re.sub("str\. ", "", description)
        description = re.sub(" complete genome sequence\.", "", description)
        description = re.sub(" complete genome\.", "", description)
        description = re.sub(" chromosome", "", description)
        description = re.sub(" DNA", "S.", description)
        description = re.sub("Merged record from ", "", description)
        description = re.sub(", wgs", "", description)
        description = re.sub("Candidatus ", "", description)
        description = re.sub(".contig.0_1, whole genome shotgun sequence.", "", description)
        description = re.sub("Protochlamydia", "P.", description)
        description = re.sub("Chlamydia", "C.", description)
        description = re.sub("Chlamydophila", "E.", description)
        description = re.sub("Estrella", "E.", description)
        description = re.sub("Rhodopirellula", "R.", description)
        description = re.sub("Methylacidiphilum", "M.", description)
        description = re.sub(" phage", "", description)
        description = re.sub("Parachlamydia", "P.", description)
        description = re.sub("Neochlamydia", "Neo.", description)
        description = re.sub("Simkania", "S.", description)
        description = re.sub("Waddlia", "W.", description)
        description = re.sub("Pirellula", "P.", description)
        description = re.sub("Rhabdochlamydiaceae sp.", "Rhabdo", description)

        taxon_id2genome_description[accession] = description
        #taxon_id2genome_description[description] = accession[0]


    with open('genome_identity_dico.json', 'r') as f:
        genome_identity = json.load(f)
    header = '\t'
    genome_identity['316'] = {}
    for i in genome_identity:
        header += '%s\t' % taxon_id2genome_description[str(i)]

    print header
    for i in genome_identity:
        row_string = '%s\t' % taxon_id2genome_description[str(i)]
        for y in genome_identity:
            #print y,i
            try:
                row_string += '%s\t' % round(genome_identity[i][y][0],2)

            except KeyError:
                try:
                    row_string += '%s\t' % round(genome_identity[y][i][0],2)
                except KeyError:
                    row_string += '-\t'
        print row_string
'''
