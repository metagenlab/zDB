#!/usr/bin/python

import manipulate_biosqldb
import gbk2circos

def locus_tag2orthogroup_size(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)
    sql='select locus_tag, count(*) from biosqldb.orthology_detail_%s group by orthogroup'

    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

def superkingdom_table():
    pass





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
        print record.id
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

def locus_tag2n_nr_hits(db_name, genome_accession):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select  A.locus_tag, B.n_hits from (select * from biosqldb.orthology_detail_%s as t1 where t1.accession="%s") A ' \
          ' left join (select t2.nr_hit_id,t2.locus_tag,t2.subject_kingdom,t2.subject_accession,count(t2.locus_tag) as n_hits ' \
          ' from blastnr.blastnr_hits_%s_%s as t2 ' \
          ' group by locus_tag) B on A.locus_tag=B.locus_tag;' %(db_name, genome_accession, db_name, genome_accession)
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
    table = pandas.DataFrame(data, columns=['locus_tag', 'protein_length','superkingdom'])
    return table



def locus_tag2n_blast_superkingdom(db_name, accession, superkingdom="Bacteria"):
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

    server, db = manipulate_biosqldb.load_db(db_name)
    data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    return data

def locus_tag2n_blast_bacterial_phylum(db_name, accession, phylum="Chlamydiae", reverse=False):
    if not reverse:
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
              ' where t4.superkingdom="Bacteria" and t4.phylum !="%s" group by locus_tag;' % (db_name,
                                                                        accession,
                                                                        db_name,
                                                                        accession,                                                                                                              db_name,
                                                                        accession,
                                                                        phylum)

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


if __name__ == '__main__':
    #print taxonomical_form('chlamydia_03_15')
    #locus_tag2best_hit_n_taxon_ids('chlamydia_03_15', 'Rhab')
    #clean_multispecies_blastnr_record('chlamydia_03_15', create_new_sql_tables=False)

    #locus_tag2n_nr_hits('chlamydia_03_15', 'Rhab')
    locus_tag2n_blast_bacteria('chlamydia_03_15', 'Rhab')