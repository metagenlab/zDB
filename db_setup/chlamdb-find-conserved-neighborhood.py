#!/usr/bin/env python



def create_locus_link_table(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create database if not exists interactions'
    server.adaptor.execute(sql,)

    sql = 'create table interactions.colocalization_table_locus_%s (locus_1 varchar(200), locus_2 varchar(200), n_links INT,' \
          'n_comparisons INT, ratio float, index locus_1 (locus_1), index locus_2 (locus_2), index n_links (n_links),' \
          ' index n_comparisons (n_comparisons), index ratio (ratio))' % biodb

    server.adaptor.execute(sql,)


def create_taxon_link_table_locus(biodb):
    # locus_a, linked_locus, ref_ortho, linked_group, taxon_a,t
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table interactions.colocalization_taxons_table_locus_%s (locus_1  varchar(200), ' \
          ' locus_2 varchar(200), group_1 varchar(200), group_2 varchar(200), taxon_1 INT,' \
          'taxon_2 INT)' % biodb

    server.adaptor.execute(sql,)




def insert_locus_links_into_mysql(biodb, locus2links):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    pair_done = []
    for locus_a in locus2links:
        for locus_b in locus2links[locus_a]:
            sql = 'insert into interactions.colocalization_table_locus_%s (locus_1, locus_2, n_links, n_comparisons, ratio)' \
                  ' values ("%s","%s",%s,%s,%s)' % (biodb,
                                                    locus_a,
                                                    locus_b,
                                                    locus2links[locus_a][locus_b][0],
                                                    locus2links[locus_a][locus_b][1],
                                                    locus2links[locus_a][locus_b][0]/float(locus2links[locus_a][locus_b][1]))
            server.adaptor.execute_and_fetchall(sql,)
    server.adaptor.commit()

def insert_taxon_links_into_mysql_locus(biodb, locus2taxons):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    for locus_1 in locus2taxons:
        for locus_2 in locus2taxons[locus_1]:
            for one_taxon_pair in locus2taxons[locus_1][locus_2]:
                    sql = 'insert into interactions.colocalization_taxons_table_locus_%s (locus_1, locus_2, group_1, ' \
                          ' group_2, taxon_1, taxon_2)' \
                          ' values ("%s","%s","%s","%s", %s, %s)' % (biodb,
                                                            locus_1,
                                                            locus_2,
                                                            one_taxon_pair[0],
                                                            one_taxon_pair[1],
                                                            one_taxon_pair[2],
                                                            one_taxon_pair[3])
                    server.adaptor.execute_and_fetchall(sql,)
    server.adaptor.commit()





def create_group_link_table(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table interactions.colocalization_table_%s (group_1 varchar(200), group_2 varchar(200), n_links INT,' \
          'n_comparisons INT, ratio float)' % biodb

    server.adaptor.execute(sql,)


def create_taxon_link_table(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table interactions.colocalization_taxons_table_%s (group_1 varchar(200), group_2 varchar(200), taxon_1 INT,' \
          'taxon_2 INT)' % biodb

    server.adaptor.execute(sql,)

def insert_group_links_into_mysql(biodb,
                      group2links):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    pair_done = []
    for i in group2links:
        for x in group2links[i]:
            if [i, x] in pair_done:
                continue
            else:
                pair_done.append([i, x])
                pair_done.append([x, i])
                sql = 'insert into interactions.colocalization_table_%s (group_1, group_2, n_links, n_comparisons, ratio)' \
                      ' values ("%s","%s",%s,%s,%s)' % (biodb,
                                                        i,
                                                        x,
                                                        group2links[i][x][0],
                                                        group2links[i][x][1],
                                                        group2links[i][x][0]/float(group2links[i][x][1]))
                server.adaptor.execute_and_fetchall(sql,)
    server.adaptor.commit()

def insert_taxon_links_into_mysql(biodb,
                      group2taxons):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    for i in group2taxons:
        for x in group2taxons[i]:
            for one_taxon_pair in group2taxons[i][x]:
                    sql = 'insert into interactions.colocalization_taxons_table_%s (group_1, group_2, taxon_1, taxon_2)' \
                          ' values ("%s","%s",%s,%s)' % (biodb,
                                                            i,
                                                            x,
                                                            one_taxon_pair[0],
                                                            one_taxon_pair[1])
                    server.adaptor.execute_and_fetchall(sql,)
    server.adaptor.commit()


def find_clusters_of_orthogroups(db_name, identity_cutoff, distance_cutoff=10000):

    '''
    ATTENTION: tous les paralogues pris en comptes
    si on a 1 prot dans genome A et 3 dans le genome B
    prot 1A sera comparee a son best hit B
    prot 1B, 2B et 2C seront comparee a 1A

    dans tous les cas cette approche est redontante car on compare tjs A vs B et B vs A...

    :param db_name: biodatabase name
    :param identity_cutoff: average ortholog identity cutoff: if genomes are too close, do not identify clusters (too much clusters)
    :param distance_cutoff: size of the considered window
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import mysqldb_plot_genomic_feature
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import pylibmc

    # for memory storage of all biorecords
    mc = pylibmc.Client(["127.0.0.1"], binary=True,
                    behaviors={"tcp_nodelay": True,
                               "ketama": True})

    server, db = manipulate_biosqldb.load_db(db_name)

    sql_locus = 'select seqfeature_id from orthology_detail_%s' % db_name
    all_locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql_locus,)]

    sql = 'select seqfeature_id, orthogroup from orthology_detail_%s' % db_name
    locus2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag,seqfeature_id from orthology_detail_%s' % db_name
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select seqfeature_id, start, stop from orthology_detail_%s' % db_name

    locus2start_end = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    orthogroup2locus_list = {}
    for locus in locus2orthogroup:
        if locus2orthogroup[locus] not in orthogroup2locus_list:
            orthogroup2locus_list[locus2orthogroup[locus]] = [locus]
        else:
            orthogroup2locus_list[locus2orthogroup[locus]].append(locus)

    sql_identity = 'select taxon_1, taxon_2, median_identity from comparative_tables.shared_og_av_id_%s' % db_name
    taxon2taxon_median_id = {}

    for row in server.adaptor.execute_and_fetchall(sql_identity,):
        if row[0] not in taxon2taxon_median_id:
            taxon2taxon_median_id[row[0]] = {}
            taxon2taxon_median_id[row[0]][row[1]] = row[2]
        else:
            taxon2taxon_median_id[row[0]][row[1]] = row[2]
    sql = 'select seqfeature_id, accession from orthology_detail_%s' % db_name
    locus2accession = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    accession_list = set(locus2accession.values())
    accession2record = {}

    locus2closest_locus_list = {}
    sql = 'select locus_1,locus_2 from comparative_tables.identity_closest_homolog_%s' % db_name
    data = server.adaptor.execute_and_fetchall(sql,)
    for i in data:
        if str(i[0]) not in locus2closest_locus_list:
            locus2closest_locus_list[str(i[0])] = [str(i[1])]
        else:
            locus2closest_locus_list[str(i[0])].append(str(i[1]))
        if i[1] not in locus2closest_locus_list:
            locus2closest_locus_list[str(i[1])] = [str(i[0])]
        else:
            locus2closest_locus_list[i[1]].append(i[0])
    # storage of all records into memory
    for accession in accession_list:
        #print accession
        rec_raw = db.lookup(accession=accession)
        try:
            new_record_reformat = mc[db_name + "_" + accession]
        except KeyError:
            #print accession, 'not in memory'
            new_record_reformat = SeqRecord(Seq(rec_raw.seq.data, rec_raw.seq.alphabet),
                                                             id=rec_raw.id, name=rec_raw.name,
                                                             description=rec_raw.description,
                                                             dbxrefs =rec_raw.dbxrefs,
                                                             features=rec_raw.features,
                                                             annotations=rec_raw.annotations)


            mc[db_name + "_" + accession]= new_record_reformat
        accession2record[accession] = new_record_reformat
    accession2taxon = manipulate_biosqldb.accession2taxon_id(server, db_name)


    # iter all orthogroups
    group2linked_groups = {}
    group2linked_taxons = {}
    all_pairs = []
    for t, ref_ortho in enumerate(list(set(locus2orthogroup.values()))):#: #: enumerate(["group_53"])
        #print t, len(list(set(locus2orthogroup.values())))
        comp_count = 0
        #reference_grp = locus2orthogroup[locus]
        group2linked_groups[ref_ortho] = {}
        locus_list = [str(i) for i in orthogroup2locus_list[ref_ortho]]

        # if no homologs, skip
        if len(locus_list) == 1:
            continue

        # iter all locus of the orthogroup
        for x, locus_a in enumerate(locus_list):
            # extract region
            print(locus2start_end)
            start_a = locus2start_end[locus_a][0]
            end_a = locus2start_end[locus_a][1]
            record = accession2record[locus2accession[locus_a]]
            size = distance_cutoff/2

            region_a = mysqldb_plot_genomic_feature.get_feature_neighborhood(start_a, end_a, record, size, 'rec')
            grp_list = []
            for feature in region_a.features:
                if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers and 'pseudogene' not in feature.qualifiers and 'translation' in feature.qualifiers:
                    locus_b = locus_tag2seqfeature_id[feature.qualifiers['locus_tag'][0]]
                    orthogroup_locus = locus2orthogroup[locus_b]
                    if orthogroup_locus not in grp_list:
                        grp_list.append(orthogroup_locus)

            # compare neighbours of all other locus to the reference locus
            try:
                closet_locus = locus2closest_locus_list[locus_a]
            except:
                continue
            for locus_b in locus_list[x+1:len(locus_list)]:

                # only consider "best hit", locus with the highest identity
                # si on a une relation 1 vs 3
                # on va avoir une seule comparaison pour le genome A vs B mais 3 pour la comparsion B vs A...
                # oubien: si plusieurs pairs: ne comparer que la paire la plus proche.
                # dans tous les cas ca va associer les groupes qui incluent de multiples paralogues.
                # keep all comparisons in memory and do it only once?
                # cas des multiples paralogues side by side
                if locus_b not in closet_locus:
                    continue

                taxon_a = accession2taxon[locus2accession[locus_a]]
                taxon_b = accession2taxon[locus2accession[locus_b]]

                # if both locus are encoded by the same taxon, skip the comparison
                if taxon_a == taxon_b:
                    continue
                try:
                    identity = taxon2taxon_median_id[taxon_a][taxon_b]
                except KeyError:
                    identity = taxon2taxon_median_id[taxon_b][taxon_a]

                # if the 2 considered genomes are too closely related, skip the comparison
                print(identity, identity_cutoff, identity < identity_cutoff)
                if identity < identity_cutoff:
                    comp_count+=1
                    #print comp_count, "comp_count"
                    start_b = locus2start_end[locus_b][0]
                    if start_b < 0:
                        start_b = 0
                    end_b = locus2start_end[locus_b][1]
                    record_b = accession2record[locus2accession[locus_b]]
                    region_b = mysqldb_plot_genomic_feature.get_feature_neighborhood(start_b, end_b, record_b, size, 'rec')
                    grp_list_b = []
                    for feature in region_b.features:
                        if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers:
                            #print feature
                            locus_b = locus_tag2seqfeature_id[feature.qualifiers['locus_tag'][0]]
                            orthogroup_locus = locus2orthogroup[locus_b]
                            if orthogroup_locus not in grp_list_b:
                                            grp_list_b.append(orthogroup_locus)


                    common = list(set(grp_list).intersection(set(grp_list_b)))
                    # remove ref group
                    try:
                        common.pop(common.index(ref_ortho))
                    except:
                        pass
                    if len(common)>0:
                        # store groups and taxons linked
                        for linked_group in common:
                            # store reciprocal relationship between the 2 genomes and the 2 groups
                            # group a vs group b == group b vs group a
                            try:
                                if [taxon_a, taxon_b] not in group2linked_taxons[ref_ortho][linked_group]:
                                    group2linked_taxons[ref_ortho][linked_group].append([taxon_a, taxon_b])
                            except KeyError:
                                try:
                                    # remove potential redundant pairs due to paralogs
                                    # all paralogs are taken into acount
                                    if [taxon_a, taxon_b] not in group2linked_taxons[linked_group][ref_ortho]:
                                        group2linked_taxons[linked_group][ref_ortho].append([taxon_a, taxon_b])
                                except KeyError:
                                    if ref_ortho in group2linked_taxons:
                                        group2linked_taxons[ref_ortho][linked_group] = [[taxon_a, taxon_b]]
                                    elif linked_group in group2linked_taxons:
                                        group2linked_taxons[linked_group][ref_ortho] = [[taxon_a, taxon_b]]
                                    else:
                                        group2linked_taxons[ref_ortho] = {}
                                        group2linked_taxons[ref_ortho][linked_group] = [[taxon_a, taxon_b]]
                            # store counts of links out of the total number of comparsions
                            if linked_group in group2linked_groups[ref_ortho]:
                                group2linked_groups[ref_ortho][linked_group][0] += 1
                            else:
                                group2linked_groups[ref_ortho][linked_group] = [1]
                    else:
                        pass
                        #print 'no common groups'
                    # check if multiple common elements
        for linked_group in group2linked_groups[ref_ortho]:
            group2linked_groups[ref_ortho][linked_group].append(comp_count)
    return group2linked_groups, group2linked_taxons


def find_clusters_of_locus(db_name, identity_cutoff, distance_cutoff=20000):

    '''
    ATTENTION: tous les paralogues pris en comptes
    si on a 1 prot dans genome A et 3 dans le genome B
    prot 1A sera comparee a son best hit B
    prot 1B, 2B et 2C seront comparee a 1A

    dans tous les cas cette approche est redontante car on compare tjs A vs B et B vs A...

    :param db_name: biodatabase name
    :param identity_cutoff: average ortholog identity cutoff: if genomes are too close, do not identify clusters (too much clusters)
    :param distance_cutoff: size of the considered window
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import mysqldb_plot_genomic_feature
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import pylibmc

    # for memory storage of all biorecords
    mc = pylibmc.Client(["127.0.0.1"], binary=True,
                    behaviors={"tcp_nodelay": True,
                               "ketama": True})

    server, db = manipulate_biosqldb.load_db(db_name)

    sql_locus = 'select seqfeature_id from orthology_detail_%s' % db_name
    all_locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql_locus,)]

    sql = 'select seqfeature_id, orthogroup from orthology_detail_%s' % db_name
    locus2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select taxon_id,description from biosqldb.bioentry;'
    taxon_id2genome_description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag,seqfeature_id from orthology_detail_%s' % db_name
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select seqfeature_id, start, stop from orthology_detail_%s' % db_name

    locus2start_end = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    orthogroup2locus_list = {}
    for locus in locus2orthogroup:
        if locus2orthogroup[locus] not in orthogroup2locus_list:
            orthogroup2locus_list[locus2orthogroup[locus]] = [locus]
        else:
            orthogroup2locus_list[locus2orthogroup[locus]].append(locus)

    sql_identity = 'select taxon_1, taxon_2, median_identity from comparative_tables.shared_og_av_id_%s' % db_name
    taxon2taxon_median_id = {}

    for row in server.adaptor.execute_and_fetchall(sql_identity,):
        if row[0] not in taxon2taxon_median_id:
            taxon2taxon_median_id[row[0]] = {}
            taxon2taxon_median_id[row[0]][row[1]] = row[2]
        else:
            taxon2taxon_median_id[row[0]][row[1]] = row[2]
    sql = 'select seqfeature_id, accession from orthology_detail_%s' % db_name
    locus2accession = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    accession_list = set(locus2accession.values())
    accession2record = {}

    locus2closest_locus_list = {}
    sql = 'select locus_1, locus_2 from comparative_tables.identity_closest_homolog2_%s' % db_name
    data = server.adaptor.execute_and_fetchall(sql,)
    for i in data:
        if str(i[0]) not in locus2closest_locus_list:
            locus2closest_locus_list[str(i[0])] = [str(i[1])]
        else:
            locus2closest_locus_list[str(i[0])].append(str(i[1]))
        #if i[1] not in locus2closest_locus_list:
        #    locus2closest_locus_list[i[1]] = [i[0]]
        #else:
        #    locus2closest_locus_list[i[1]].append(i[0])
    # storage of all records into memory
    for accession in accession_list:
        #print accession
        rec_raw = db.lookup(accession=accession)
        try:
            new_record_reformat = mc[db_name + "_" + accession]
            print(accession, 'in memory')
        except KeyError:
            print(accession, 'NOT in memory')
            new_record_reformat = SeqRecord(Seq(rec_raw.seq.data, rec_raw.seq.alphabet),
                                                             id=rec_raw.id,
                                                             name=rec_raw.name,
                                                             description=rec_raw.description,
                                                             dbxrefs =rec_raw.dbxrefs,
                                                             features=rec_raw.features,
                                                             annotations=rec_raw.annotations)


            mc[db_name + "_" + accession]= new_record_reformat
        accession2record[accession] = new_record_reformat
    accession2taxon = manipulate_biosqldb.accession2taxon_id(server, db_name)

    # iter all orthogroups
    locus2linked_locus = {}

    locus2linked_taxons = {}
    all_pairs = []

    # iterate all orthogroups
    for t, ref_ortho in enumerate(list(set(locus2orthogroup.values()))):#: #: enumerate(["group_53"])
        print ('group %s / %s' % (t, len(list(set(locus2orthogroup.values())))))

        tmp_dico = {}

        # locus list of the considered group
        locus_list = orthogroup2locus_list[ref_ortho]

        # if a single locus, slip
        if len(locus_list) == 1:
            continue

        # iter all locus of the orthogroup
        for x, locus_a in enumerate(locus_list):
            locus_a = str(locus_a)
            locus2linked_locus[locus_a] = {}
            tmp_dico[locus_a] = {}
            locus2linked_taxons[locus_a] = {}
            # for each locus, initiate the count of comparisons
            comp_count = 0

            # extract region
            start_a = locus2start_end[locus_a][0]
            end_a = locus2start_end[locus_a][1]
            record = accession2record[locus2accession[locus_a]]
            size = distance_cutoff/2

            region_a = mysqldb_plot_genomic_feature.get_feature_neighborhood(start_a, end_a, record, size, 'rec')

            # get list of orthogroups in the neiborhood & get corresp between locus and groups
            grp_list = []
            grp2locus = {}
            for feature in region_a.features:
                if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers:
                    locus_b = str(locus_tag2seqfeature_id[feature.qualifiers['locus_tag'][0]])
                    orthogroup_locus = locus2orthogroup[locus_b]
                    if orthogroup_locus not in grp2locus:
                        grp2locus[orthogroup_locus] = [locus_b]
                    else:
                        grp2locus[orthogroup_locus].append(locus_b)
                    if orthogroup_locus not in grp_list:
                        grp_list.append(orthogroup_locus)

            # compare neighbours of all other locus to the reference locus
            try:
                closet_locus = locus2closest_locus_list[locus_a]
            except KeyError:
                continue
            for locus_b in closet_locus:#locus_list[x+1:len(locus_list)]:

                taxon_a = accession2taxon[locus2accession[locus_a]]
                taxon_b = accession2taxon[locus2accession[locus_b]]

                # if both locus are encoded by the same taxon, skip the comparison
                if taxon_a == taxon_b:
                    continue

                try:
                    identity = taxon2taxon_median_id[taxon_a][taxon_b]
                except KeyError:
                    identity = taxon2taxon_median_id[taxon_b][taxon_a]

                # if the 2 considered genomes are too closely related, skip the comparison
                if identity < identity_cutoff:
                    # increment the number of effective comparisons
                    comp_count += 1
                    start_b = locus2start_end[locus_b][0]
                    # if border of contig/chromosome
                    if start_b < 0:
                        start_b = 0
                    end_b = locus2start_end[locus_b][1]
                    record_b = accession2record[locus2accession[locus_b]]
                    region_b = mysqldb_plot_genomic_feature.get_feature_neighborhood(start_b, end_b, record_b, size, 'rec')

                    # get group list b
                    grp_list_b = []
                    for feature in region_b.features:
                        if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers:
                            locus_b = str(locus_tag2seqfeature_id[feature.qualifiers['locus_tag'][0]])
                            orthogroup_locus = locus2orthogroup[locus_b]
                            if orthogroup_locus not in grp_list_b:
                                            grp_list_b.append(orthogroup_locus)

                    # get list of common groups
                    common = list(set(grp_list).intersection(set(grp_list_b)))
                    # remove ref group
                    try:
                        common.pop(common.index(ref_ortho))
                    except:
                        with open('problems.txt', 'a') as f:
                            f.write('%s\t%s\t%s\n' % (ref_ortho, locus_a,locus_b))
                    if len(common)>0:
                        # store locus and taxons linked
                        for linked_group in common:
                            # store reciprocal relationship between the 2 genomes and the 2 locus
                            # we can have more than one locus/group (i.e identical genes side by side)
                            for linked_locus in grp2locus[linked_group]:
                                # if reverse comparison was already made
                                #if linked_locus in locus2linked_locus:
                                #    continue
                                if linked_locus not in tmp_dico[locus_a]: #locus2linked_locus
                                    # store n link and n comparisons
                                    #locus2linked_locus[locus_a][linked_locus] = [1, '-']
                                    tmp_dico[locus_a][linked_locus] = [1, '-']
                                    locus2linked_taxons[locus_a][linked_locus] = [[ref_ortho, linked_group, taxon_a,taxon_b]]
                                else:
                                    #locus2linked_locus[locus_a][linked_locus][0] += 1
                                    tmp_dico[locus_a][linked_locus][0] += 1
                                    locus2linked_taxons[locus_a][linked_locus].append([ref_ortho, linked_group, taxon_a,taxon_b])

            # end of loop for locus_a: store the number of comparisons done
            for linked_locus in tmp_dico[locus_a]: # locus2linked_locus
                #locus2linked_locus[locus_a][linked_locus][1] = comp_count
                tmp_dico[locus_a][linked_locus][1] = comp_count
            #if len(locus2linked_locus[locus_a]) == 0:
            #    del locus2linked_locus[locus_a]
            if len(tmp_dico[locus_a]) > 0:
                #print 'insert!'
                #print tmp_dico
                for locus_b in tmp_dico[locus_a]:
                    # only add minimum of 50% links
                    if tmp_dico[locus_a][locus_b][0]/float(tmp_dico[locus_a][locus_b][1]) > 0.5:
                        sql = 'insert into interactions.colocalization_table_locus_%s (locus_1, locus_2, n_links, n_comparisons, ratio)' \
                              ' values ("%s","%s",%s,%s,%s)' % (db_name,
                                                                locus_a,
                                                                locus_b,
                                                                tmp_dico[locus_a][locus_b][0],
                                                                tmp_dico[locus_a][locus_b][1],
                                                                tmp_dico[locus_a][locus_b][0]/float(tmp_dico[locus_a][locus_b][1]))
                        server.adaptor.execute(sql,)
                        server.adaptor.commit()

            if len(locus2linked_taxons[locus_a]) == 0:
                del locus2linked_taxons[locus_a]
            #print locus2linked_locus
    return locus2linked_locus, locus2linked_taxons

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--database_name', type=str, help="Database name")

    args = parser.parse_args()

    #group_links, taxon_links = find_clusters_of_orthogroups('chlamydia_04_16', 60, 10000)
    #create_group_link_table('chlamydia_04_16')
    #create_taxon_link_table('chlamydia_04_16')
    #insert_group_links_into_mysql('chlamydia_04_16', group_links)
    #insert_taxon_links_into_mysql('chlamydia_04_16', taxon_links)

    create_locus_link_table(args.database_name)
    create_taxon_link_table_locus(args.database_name)
    locus_links, taxon_links = find_clusters_of_locus(args.database_name, 60, 20000)

    insert_locus_links_into_mysql(args.database_name, locus_links)
    insert_taxon_links_into_mysql_locus(args.database_name, taxon_links)
