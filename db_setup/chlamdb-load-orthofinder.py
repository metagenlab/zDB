#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-

from chlamdb.orthomcl.parse_mclOtput import parse_orthomcl_output
from chlamdb.biosqldb import manipulate_biosqldb
import os
import numpy as np
import numpy
from multiprocessing import Process, Queue
import time
from multiprocessing import cpu_count
import MySQLdb


def add_orthogroup_term(server):
    # => ajouter orthogroup � la liste des term_id si n'existe pas
    # | term_id | name                   | definition | identifier | is_obsolete | ontology_id |
    # |      21 | CDS                    | NULL       | NULL       | NULL        |           2 |

    # seqfeature_qualifier_value
    # | term_id | seqfeature_id | rank | value | name | definition | identifier | is_obsolete | ontology_id |
    # |      16 |          5498 |    1 | hoxN  | gene | NULL       | NULL       | NULL        |           1 |
    try:
        # test if orthogroup term exist
        sql1 = 'SELECT term_id FROM term WHERE NAME = "orthogroup";'
        return server.adaptor.execute_and_fetchall(sql1)[0][0]
    except:
        #print "adding orthogroup seqfeature term to", server
        # add orthogroup, ontology term: 2
        sql1 = 'select ontology_id from ontology where name="SeqFeature Keys"'
        id = server.adaptor.execute_and_fetchall(sql1)[0][0]
        sql2 = 'INSERT INTO term (name, ontology_id) values ("orthogroup", %s);' % id
        server.adaptor.execute(sql2)
        server.adaptor.commit()
    return server.adaptor.execute_and_fetchall(sql1)[0][0]


def add_orthogroup_to_seq(server, 
                          locus_tag2orthogroup, 
                          locus_tag2seqfeature_id_dico):
    #| seqfeature_id | term_id | rank | value   |
    rank = 1
    term_id = add_orthogroup_term(server)
    for n, locus_tag in enumerate(locus_tag2orthogroup.keys()):
        if n % 100 == 0:
            print("%s / %s" % (n, len(locus_tag2orthogroup.keys())))
        seqfeature_id = locus_tag2seqfeature_id_dico[locus_tag]
        group = locus_tag2orthogroup[locus_tag]
        sql = 'INSERT INTO seqfeature_qualifier_value (seqfeature_id, term_id, `rank`, value) values (%s, %s, %s, "%s");' % (seqfeature_id,
                                                                                                                           term_id,
                                                                                                                           rank,
                                                                                                                           group)
        server.adaptor.execute(sql)
    server.adaptor.commit()


def create_indexes(server):

    # orthology_orthogroup    
    sql_index1 = 'create index ooo on orthology_orthogroup(orthogroup_name)'
    
    # orthology_seqfeature_id2orthogroup
    sql_index2 = 'create index osoi on orthology_seqfeature_id2orthogroup(seqfeature_id)'
    sql_index3 = 'create index ooos on orthology_seqfeature_id2orthogroup(orthogroup_id)'
    
    # annotation_seqfeature_id2locus
    sql_index4 = 'create index asls on annotation_seqfeature_id2locus(seqfeature_id)'
    sql_index5 = 'create index aslb on annotation_seqfeature_id2locus(bioentry_id)'
    sql_index6 = 'create index asll on annotation_seqfeature_id2locus(locus_tag)'
    sql_index7 = 'create index aslt on annotation_seqfeature_id2locus(taxon_id)'
    
    # annotation_seqfeature_id2CDS_annotation
    sql_index8 = 'create index asics on annotation_seqfeature_id2CDS_annotation(seqfeature_id)'
    
    # annotation_seqfeature_id2RNA_annotation
    sql_index9 = 'create index asirs on annotation_seqfeature_id2RNA_annotation(seqfeature_id)'
    
    # orthology_detail
    sql_index10 = 'create index odo on orthology_detail(orthogroup)'
    sql_index11 = 'create index odt on orthology_detail(taxon_id)'
    sql_index12 = 'create index odl on orthology_detail(locus_tag)'
    sql_index13 = 'create index ods on orthology_detail(seqfeature_id)'
    
    
    server.adaptor.execute(sql_index1)
    server.adaptor.execute(sql_index2)
    server.adaptor.execute(sql_index3)
    server.adaptor.execute(sql_index4)
    server.adaptor.execute(sql_index5)
    server.adaptor.execute(sql_index6)
    server.adaptor.execute(sql_index7)
    server.adaptor.execute(sql_index8)
    server.adaptor.execute(sql_index9)
    server.adaptor.execute(sql_index10)
    server.adaptor.execute(sql_index11)
    server.adaptor.execute(sql_index12)
    server.adaptor.execute(sql_index13)
    
    
    
def create_orthology_tables(server):
    sql = 'CREATE TABLE if not exists orthology_orthogroup (orthogroup_id INTEGER PRIMARY KEY,' \
          ' orthogroup_name varchar(400),' \
          ' orthogroup_size INT,' \
          ' n_genomes INT)'
          
    server.adaptor.execute(sql)

    sql0 = 'CREATE TABLE if not exists orthology_seqfeature_id2orthogroup (seqfeature_id INT,' \
           ' orthogroup_id INT)'

    server.adaptor.execute(sql0)

    sql1 = 'CREATE TABLE if not exists annotation_seqfeature_id2locus (seqfeature_id INT,' \
          ' feature_type_id INT,' \
          ' taxon_id INT,' \
          ' pseudogene INT, ' \
          ' bioentry_id VARCHAR(100) NOT NULL, ' \
          ' locus_tag VARCHAR(100) NOT NULL, ' \
          ' start INT, ' \
          ' stop INT, ' \
          ' strand INT)'
          
    server.adaptor.execute(sql1)

    sql2 = 'CREATE TABLE if not exists annotation_seqfeature_id2CDS_annotation (seqfeature_id INT,' \
          ' gene VARCHAR(100) NOT NULL, ' \
          ' protein_id VARCHAR(100) NOT NULL, ' \
          ' product TEXT NOT NULL, ' \
          ' translation TEXT NOT NULL, ' \
          ' SP INT,' \
          ' TM INT)'

    server.adaptor.execute(sql2)

    sql3 = 'CREATE TABLE if not exists annotation_seqfeature_id2RNA_annotation (seqfeature_id INT,' \
          ' product TEXT NOT NULL)'

    server.adaptor.execute(sql3)    


def create_annotation_seqfeature_id2locus_and_cds(server,
                                                  locus_tag2seqfeature_id_dico,
                                                  locus_tag2taxon_dico,
                                                  seqfeature_id2protein_id_dico,
                                                  seqfeature_id2gene_dico,
                                                  seqfeature_id2product_dico,
                                                  seqfeature_id2translation_dico,
                                                  seqfeature_id2bioentry_id_dico,
                                                  seqfeature2location_dico,
                                                  seqfeature_id2feature_type_id,
                                                  pseudogene_feature_list):

    # dico locus type (CDS, rRNA, tRNA)
    sql = 'select term_id, name from term'
    term_id2term = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # fill seqfeature_id2CDS_annotation_ and seqfeature_id2locus_
    #print 'Fill CDS tables (CDS and RNA tables)'
    #print 'pseudo list:', pseudogene_feature_list[0:10]
    #print 'seqfeature_id2bioentry_id_dico', seqfeature_id2bioentry_id_dico.keys()[0:10]
    #print 'seqfeature_id2feature_type_id', seqfeature_id2feature_type_id.keys()[0:10]
    #print 'seqfeature_id2protein_id_dico', seqfeature_id2protein_id_dico.keys()[0:10]

    for n, locus_tag in enumerate(locus_tag2seqfeature_id_dico):
        if n % 5000 == 0:
            print("annotation_seqfeature_id2locus: %s / %s" % (n, len(locus_tag2seqfeature_id_dico)))
        seqfeature_id = locus_tag2seqfeature_id_dico[locus_tag]
        taxon_id = locus_tag2taxon_dico[locus_tag]
        start = seqfeature2location_dico[seqfeature_id][0]
        end = seqfeature2location_dico[seqfeature_id][1]
        strand = seqfeature2location_dico[seqfeature_id][2]
        bioentry_id = seqfeature_id2bioentry_id_dico[str(seqfeature_id)]
        seqfeature_type_id = seqfeature_id2feature_type_id[seqfeature_id]

        if seqfeature_id not in pseudogene_feature_list:
            if term_id2term[str(seqfeature_type_id)] == 'CDS':
                try:
                    translation = seqfeature_id2translation_dico[str(seqfeature_id)]
                    pseudo = 0
                except KeyError:
                    print("Missing translation for locus: %s, consider it as pseudogene" % locus_tag)
                    pseudo = 1
        else:
            pseudo = 1

        sql = 'insert into annotation_seqfeature_id2locus (seqfeature_id, feature_type_id, taxon_id, ' \
              ' pseudogene, bioentry_id, locus_tag, start, stop, strand) ' \
              ' values (%s, %s, %s, %s, %s, "%s", %s, %s, %s)' % (seqfeature_id,
                                                                  seqfeature_type_id,
                                                                  taxon_id,
                                                                  pseudo,
                                                                  bioentry_id,
                                                                  locus_tag,
                                                                  start,
                                                                  end,
                                                                  strand)


        server.adaptor.execute(sql,)

        # include rRNA, tRNA
        if pseudo != 1:
            if term_id2term[str(seqfeature_type_id)] == 'CDS':
                try:
                    protein_id = seqfeature_id2protein_id_dico[str(seqfeature_id)]
                except KeyError:
                    protein_id = locus_tag
                try:
                    gene = seqfeature_id2gene_dico[str(seqfeature_id)]
                except KeyError:
                    gene = "-"

                try:
                    product = seqfeature_id2product_dico[str(seqfeature_id)]
                except KeyError:
                    product = "-"

                sql = 'insert into annotation_seqfeature_id2CDS_annotation (seqfeature_id, ' \
                      ' gene, protein_id, product, translation, SP, TM) ' \
                      ' values (%s, "%s", "%s", "%s", "%s", %s, %s)' % (seqfeature_id,
                                                                        gene,
                                                                        protein_id,
                                                                        product,
                                                                        translation,
                                                                        0,
                                                                        0)
                ##print sql
                server.adaptor.execute(sql,)

            elif term_id2term[str(seqfeature_type_id)] in ['rRNA', 'tRNA']:
                try:
                    product = seqfeature_id2product_dico[str(seqfeature_id)]
                except KeyError:
                    product = "-"
                sql = 'insert into annotation_seqfeature_id2RNA_annotation ' \
                      ' (seqfeature_id, product)' \
                      ' values (%s, "%s")' % (seqfeature_id,
                                              product)
                server.adaptor.execute(sql,)
            else:
                print ('Unkown feature type', seqfeature_type_id)
        server.adaptor.commit()


def create_orthology_orthogroup(server,
                                group2orthogroup_size,
                                group2family_size,
                                orthomcl_groups2locus_tag_list):

    # create orthogroup table
    #print 'Fill orthogroup table'
    for i in range(0, len(orthomcl_groups2locus_tag_list.keys())):
        if i % 1000 == 0:
            print("orthology_orthogroup: %s/%s groups" % (i, len(orthomcl_groups2locus_tag_list)))
        group = "group_%s" % i

        try:
            orthogroup_size = group2orthogroup_size[group]
        except:
            #print 'no size?----------', group
            continue
        n_genomes = group2family_size[group]

        sql = 'insert into orthology_orthogroup (orthogroup_name, orthogroup_size, n_genomes) ' \
              ' values ("%s", %s, %s)' % (group,
                                          orthogroup_size,
                                          n_genomes)

        server.adaptor.execute(sql,)

    server.adaptor.commit()


    
def orthology_seqfeature_id2orthogroup(server,
                                        orthomcl_groups2locus_tag_list,
                                        locus_tag2seqfeature_id_dico):


    # get orthogroup ids
    sql = 'select orthogroup_name, orthogroup_id from orthology_orthogroup'

    orthogroup2orthogroup_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # fill seqfeature_id2orthogroup
    #print 'fill seqfeature_id2orthogroup'
    for i in range(0, len(orthomcl_groups2locus_tag_list.keys())):
        if i % 1000 == 0:
            print("orthology_seqfeature_id2orthogroup: %s/%s groups" % (i, len(orthomcl_groups2locus_tag_list)))
        group = "group_%s" % i

        locus_tag_list = orthomcl_groups2locus_tag_list[group]
        orthogroup_id = orthogroup2orthogroup_id[group]

        # iterate proteins of the group, insert in the seqfeature2orthogroup table
        for locus_tag in locus_tag_list:
            # if the protein has no protein id (and a locus tag only => prokka annotation)
            seqfeature_id = locus_tag2seqfeature_id_dico[locus_tag]

            sql = 'INSERT INTO orthology_seqfeature_id2orthogroup(seqfeature_id, orthogroup_id) ' \
                  ' values(%s,%s);' % (seqfeature_id,
                                       orthogroup_id)
            print(sql)
            server.adaptor.execute(sql)
        server.adaptor.commit()


def get_all_orthogroup_size(server, biodatabase_name):
    """
    return a dictonary with orthogroup size"
    """

    sql = ' select seqfeature_qualifier_value.value, COUNT(*) from seqfeature_qualifier_value' \
          ' inner join term on seqfeature_qualifier_value.term_id = term.term_id and term.name = "orthogroup"' \
          ' inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id' \
          ' inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id' \
          ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "%s"' \
          ' group by seqfeature_qualifier_value.value' % biodatabase_name

    result = server.adaptor.execute_and_fetchall(sql,)

    return manipulate_biosqldb._to_dict(result)


def get_family_size(server, biodatabase_name):
    '''
    return number of genomes in which the family is present (multiple copies are only counted ones)
    '''
    table = manipulate_biosqldb.get_orthology_table2(server, biodatabase_name)
    groupid2family_size = {}
    for row in table:
        #
        family = 0
        # todo really -1?????
        # row[1:] because first column is the orthogroup name
        for i in row[1:]:
            if i >0:
                family+=1
        groupid2family_size[row[0]] = family
    return groupid2family_size


def create_orthogroup_table_legacy(server, 
                                   biodatabase_name,
                                   orthomcl_groups2locus_tag_list,
                                   locus_tag2seqfeature_id_dico,
                                   protein_id2seqfeature_id_dico,
                                   locus_tag2taxon_dico,
                                   protein_id2taxon_dico,
                                   locus_tag2accession_dico,
                                   protein_id2accession_dico,
                                   seqfeature_id2locus_tag_dico,
                                   seqfeature_id2protein_id_dico,
                                   seqfeature_id2gene_dico,
                                   seqfeature_id2product_dico,
                                   seqfeature_id2translation_dico,
                                   seqfeature_id2organism_dico,
                                   seqfeature2location_dico,
                                   group2orthogroup_size,
                                   group2family_size):
    import re

    sql = 'CREATE TABLE orthology_detail (orthogroup VARCHAR(100) NOT NULL, ' \
          ' taxon_id INT, ' \
          ' accession VARCHAR(100) NOT NULL, ' \
          ' locus_tag VARCHAR(100) NOT NULL, ' \
          ' protein_id VARCHAR(100) NOT NULL, ' \
          ' start INT, ' \
          ' stop INT, ' \
          ' strand INT, ' \
          ' gene VARCHAR(100) NOT NULL, ' \
          ' product TEXT NOT NULL, ' \
          ' translation TEXT NOT NULL, ' \
          ' organism TEXT NOT NULL, ' \
          ' orthogroup_size INT,' \
          ' n_genomes INT,' \
          ' seqfeature_id INT)'

    server.adaptor.execute(sql)

    for i in range(0, len(orthomcl_groups2locus_tag_list.keys())):
        if i % 1000 == 0:
            print("%s / %s" % (i, len(orthomcl_groups2locus_tag_list)))
        group = "group_%s" % i
        try:
            orthogroup_size = group2orthogroup_size[group]
        except KeyError:
            continue
        n_genomes = group2family_size[group]
        locus_tag_list = orthomcl_groups2locus_tag_list[group]

        # iterate all proteins from orthomcl file
        # add data to the summary table othology_detail...
        for locus_tag in locus_tag_list:
            # if the protein has no protein id (and a locus tag only => prokka annotation)
            try:
                seqfeature_id = locus_tag2seqfeature_id_dico[locus_tag]
            except KeyError:
                print (locus_tag, 'not loaded')
                continue
            taxon_id = locus_tag2taxon_dico[locus_tag]
            accession = locus_tag2accession_dico[locus_tag]
            start = seqfeature2location_dico[seqfeature_id][0]
            end = seqfeature2location_dico[seqfeature_id][1]
            strand = seqfeature2location_dico[seqfeature_id][2]
            organism = seqfeature_id2organism_dico[str(seqfeature_id)]
            translation = seqfeature_id2translation_dico[str(seqfeature_id)]

            try:
                gene = seqfeature_id2gene_dico[str(seqfeature_id)]
            except KeyError:
                gene = "-"

            try:
                product = seqfeature_id2product_dico[str(seqfeature_id)]
            except KeyError:
                product = "-"

            seqfeature_id = locus_tag2seqfeature_id_dico[locus_tag]
            try:
                protein_id = seqfeature_id2protein_id_dico[str(seqfeature_id)]
            except KeyError:
                protein_id = locus_tag

            sql = 'INSERT INTO orthology_detail (orthogroup, taxon_id, accession, locus_tag, protein_id, start, ' \
                  'stop, strand, gene, product, translation, organism, orthogroup_size, n_genomes, seqfeature_id) ' \
                  ' values ("%s", %s, "%s", "%s", "%s", %s, %s, %s, "%s", "%s", "%s", "%s", %s, %s, %s);' % (group, taxon_id,
                                                                                                             accession,
                                                                                                             locus_tag,
                                                                                                             protein_id,
                                                                                                             start, end,
                                                                                                             strand,
                                                                                                             gene,
                                                                                                             product,
                                                                                                             translation,
                                                                                                             organism,
                                                                                                             orthogroup_size,
                                                                                                             n_genomes,
                                                                                                             seqfeature_id)
            try:
                server.adaptor.execute(sql)
            except:
                print('problem with:')
                print(sql)
    server.adaptor.commit()


def plot_orthogroup_size_distrib(server, biodatabase_name, out_name = "orthogroup_size_distrib.svg", taxon_id = False):

    from collections import Counter
    import numpy as np
    import pandas
    import itertools
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import math
    import re

    if not taxon_id:

        all_taxons = manipulate_biosqldb.get_taxon_id_list(server, biodatabase_name)
        for key in all_taxons:
            sql = 'select description from bioentry where taxon_id=%s and description not like "%%%%plasmid%%%%" limit 1;'
            #print sql, key
            description = server.adaptor.execute_and_fetchall(sql, key)[0][0]

        #print all_taxons
        all_data = {}
        for taxon in all_taxons:
            sql = 'select orthogroup,`%s` from comparative_tables_orthology where `%s` >0' % (taxon, 
                                                                                              taxon)
            result = manipulate_biosqldb._to_dict(server.adaptor.execute_and_fetchall(sql, ))
            all_data[taxon] = Counter(result.values())

        # get maximum group size and max number of proteins for a single group
        max_orthogroup_size = 0
        max_n_groups = 0
        for key in all_data.keys():
            dico = all_data[key]
            temp = int(max(dico.keys()))
            temp2 = int(max(dico.values()))
            if  temp > max_orthogroup_size:
                max_orthogroup_size = temp
            if temp2 > max_n_groups:
                max_n_groups = temp2

        #for key in all_data.keys():
        #    dico= all_data[key]
        #    for i in range(1, max_orthogroup_size):
        #        if i not in dico.keys():
        #            dico[i] = 0

        # transform into pandas serie objects
        complete_data = {}
        for key in all_data.keys():
            dico = all_data[key]
            df = pandas.Series(dico)
            complete_data[key] = df

        # create figure with one subplot per genome
        # plot per genome group size distribution
        n = int(math.ceil(len(all_taxons)/2.0))
        #print n
        if n%2 ==0:
            n+=1
        fig, axes = plt.subplots(nrows=n, ncols=2)
        for key, location in zip(complete_data.keys(), list(itertools.product(range(n),repeat=2))):
            serie = complete_data[key]
            sql = 'select description from bioentry where taxon_id=%s and description not like "%%%%plasmid%%%%" limit 1;'
            #print sql, key
            description = server.adaptor.execute_and_fetchall(sql, key)[0][0]
            #print description

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
            p = serie.plot(ax=axes[location[1],location[0]],kind="bar", alpha=1, color='r')
            p.set_ylim(-100, max_n_groups + 0.1*max_n_groups)
            p.set_title(description)

            for rect in p.patches:
                height= rect.get_height()
                p.text(rect.get_x()+rect.get_width()/2., height + 0.02*max_n_groups, int(height), ha='center', va='bottom')

        # plot all genomes group size distribution
        #server, db = manipulate_biosqldb.load_db(db_name)
        sql='select orthogroup, count(*) from orthology_detail group by orthogroup'

        # based on the new orthology detail table
        all_grp_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        # based on biosqldb tables
        #all_grp_size = get_all_orthogroup_size(server, biodatabase_name)

        df = pandas.Series(Counter(all_grp_size.values()))
        p = df.plot(kind="bar", alpha=0.5, color = 'b')
        p.set_ylim(-100, max(df) + 0.1*max(df))
        p.set_title("Orthogroup size distribution")

        # add group size as text at the top of the bars
        for rect in p.patches:
            height= rect.get_height()
            p.text(rect.get_x()+rect.get_width()/2., height + 0.02*max_n_groups, int(height), ha='center', va='bottom')
        fig.set_size_inches(18,6*n)
        plt.savefig(out_name, format="svg")


    if taxon_id:
        #print "ID!!!"
        sql = 'select orthogroup,`%s` from %s where `%s` >0'
        othogroup_id2size = manipulate_biosqldb._to_dict(server.adaptor.execute_and_fetchall(sql, (taxon_id, biodatabase_name, taxon_id)))

        serie = pandas.Series(Counter(othogroup_id2size.values()))
        p = serie.plot(kind="bar", alpha=0.5, color = 'b')
        p.set_ylim(-100, max(serie) + 0.1*max(serie))
        sql = 'select description from bioentry where taxon_id=%s limit 1;'
        description = server.adaptor.execute_and_fetchall(sql, taxon_id)
        p.set_title(description[0][0])

        for rect in p.patches:
            height= rect.get_height()
            p.text(rect.get_x()+rect.get_width()/2., height + 0.02*max(serie), int(height), ha='center', va='bottom')
        #plt.show()
        plt.savefig(out_name, format="svg")

def get_orthology_matrix_separate_plasmids(server, biodatabase_name):
    genomes = manipulate_biosqldb.get_genome_description_list(server, biodatabase_name)
    for genome in genomes:
        #print genome
        dico = manipulate_biosqldb.bioentry_name2orthoup_size(server, biodatabase_name, genome)
        #print dico



def get_orthology_matrix_merging_plasmids_biosqldb(server, biodatabase_name):

    all_orthogroups = get_all_orthogroup_size(server, biodatabase_name)
    #server, db = manipulate_biosqldb.load_db(biodatabase_name)
    #sql='select orthogroup, count(*) from orthology_detail group by orthogroup' % biodatabase_name

    # based on the new orthology detail table
    #all_orthogroups = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


    n_groups = len(all_orthogroups.keys())
    print ("ngroups", n_groups)


    all_taxons = manipulate_biosqldb.get_taxon_id_list(server, biodatabase_name)
    ##print all_taxons
    #print 'get dico ortho count'
    detailed_orthology_count = {}
    for group in all_orthogroups.keys():
        detailed_orthology_count[group] = {}
        for taxon_id in all_taxons:
            detailed_orthology_count[group][int(taxon_id)]= 0


    for taxon_id in all_taxons:
        ##print taxon_id
        dico = manipulate_biosqldb.taxon_id2orthogroup_size(server, biodatabase_name, taxon_id)
        ##print dico
        #import time
        #time.sleep(3)
        ##print "taxon", taxon_id
        for group in dico.keys():
            detailed_orthology_count[group][int(taxon_id)] += dico[group]
    #print 'count ok'
    ##print detailed_orthology_count
    return detailed_orthology_count



def get_orthology_matrix_merging_plasmids_own_tables(server, biodatabase_name):

    #all_orthogroups = get_all_orthogroup_size(server, biodatabase_name)
    server, db = manipulate_biosqldb.load_db(biodatabase_name)
    sql='select orthogroup, count(*) from orthology_detail group by orthogroup'

    # based on the new orthology detail table
    all_orthogroups = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


    n_groups = len(all_orthogroups.keys())
    ##print "ngroups", n_groups


    all_taxons = manipulate_biosqldb.get_taxon_id_list(server, biodatabase_name)
    ##print all_taxons
    #print 'get dico ortho count'
    detailed_orthology_count = {}
    for group in all_orthogroups.keys():
        detailed_orthology_count[group] = {}
        for taxon_id in all_taxons:
            detailed_orthology_count[group][int(taxon_id)]= 0


    for taxon_id in all_taxons:
        #print taxon_id
        dico = manipulate_biosqldb.taxon_id2orthogroup_size(server, biodatabase_name, taxon_id)
        #print dico
        #import time
        #time.sleep(3)
        #print "taxon", taxon_id
        for group in dico.keys():
            detailed_orthology_count[group][int(taxon_id)] += dico[group]
    #print 'count ok'
    #print detailed_orthology_count
    return detailed_orthology_count



def create_orthology_mysql_table(server, orthogroup2detailed_count, biodatabase_name):

    """
    CREATE TABLE table_name (column_name column_type);
    """
    #print
    taxon_id_list = list(orthogroup2detailed_count[list(orthogroup2detailed_count.keys())[0]].keys())

    sql = "CREATE TABLE comparative_tables_orthology (orthogroup VARCHAR(100) NOT NULL"
    for i in taxon_id_list:
        sql+=" ,`%s` INT" % i
    sql+=")"
    server.adaptor.execute(sql)

    for group in list(orthogroup2detailed_count.keys()):
        values='"%s", ' % group
        columns='orthogroup, '
        taxons = list(orthogroup2detailed_count[group].keys())
        for i in range(0, len(taxons)-1):
            values += " %s," % orthogroup2detailed_count[group][taxons[i]]
            columns+="`%s`, " % taxons[i]
        values += "%s" % orthogroup2detailed_count[group][taxons[-1]]
        columns += "`%s`" % taxons[-1]
        sql = "INSERT INTO comparative_tables_orthology (%s) values (%s);" %  (columns, values)
        #print sql
        server.adaptor.execute(sql)
        server.adaptor.commit()

def get_taxon_id_description_table(server, biodatabase_name, taxon_id_list):

    template = "taxon_id = %s "
    query = "where "
    for i in range(0, len(taxon_id_list)-1):
        query += template % taxon_id_list[i] + "or "
    query += template % taxon_id_list[-1]

    sql = 'select accession, taxon_id, bioentry.description from bioentry' \
' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "%s" %s;'
    sql = sql % (biodatabase_name, query)
    result = server.adaptor.execute_and_fetchall(sql, )
    return result

def _chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]


def get_non_protein_coding_locus(server, biodb,term_list=["rRNA","tRNA","assembly_gap"]):


    sql0 = 'CREATE table non_protein_coding_locus (seqfeature_id INT, ' \
           ' taxon_id INT, ' \
           ' bioentry_id INT, ' \
           ' accession varchar(200), ' \
           ' type varchar(200), ' \
           ' locus_tag varchar(200), ' \
           ' product TEXT, ' \
           ' start INT, ' \
           ' stop INT, ' \
           ' strand INT)' % biodb

    server.adaptor.execute(sql0,)

    filter = '"'+'","'.join(term_list)+'"'
    #print filter

    sql = 'select t2.bioentry_id,taxon_id,accession,seqfeature_id,type_term_id,t4.name ' \
          ' from biodatabase as t1 inner join bioentry as t2 on t1.biodatabase_id =t2.biodatabase_id ' \
          ' inner join seqfeature as t3 on t2.bioentry_id=t3.bioentry_id ' \
          ' inner join term t4 on t3.type_term_id=t4.term_id where t1.name="%s" and t4.name in (%s);' % (biodb, filter)
    #print sql

    non_protein_coding_data = server.adaptor.execute_and_fetchall(sql,)

    sql2 = 'select A.seqfeature_id,value from ' \
           ' (select t2.bioentry_id,taxon_id,accession,seqfeature_id,type_term_id,t4.name ' \
           ' from biodatabase as t1 inner join bioentry as t2 on t1.biodatabase_id =t2.biodatabase_id ' \
           ' inner join seqfeature as t3 on t2.bioentry_id=t3.bioentry_id ' \
           ' inner join term t4 on t3.type_term_id=t4.term_id where t1.name="%s" ' \
           ' and t4.name in (%s))A inner join seqfeature_qualifier_value B on A.seqfeature_id=B.seqfeature_id ' \
           'inner join term as C on B.term_id=C.term_id where C.name="locus_tag"' % (biodb, filter)
    #print sql2
    seqfeature_id2locus_tag = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    #print "seqfeature_id2locus_tag", seqfeature_id2locus_tag
    sql3 = 'select A.seqfeature_id,value from ' \
           ' (select t2.bioentry_id,taxon_id,accession,seqfeature_id,type_term_id,t4.name ' \
           ' from biodatabase as t1 inner join bioentry as t2 on t1.biodatabase_id =t2.biodatabase_id ' \
           ' inner join seqfeature as t3 on t2.bioentry_id=t3.bioentry_id ' \
           ' inner join term t4 on t3.type_term_id=t4.term_id where t1.name="%s" ' \
           ' and t4.name in (%s))A inner join seqfeature_qualifier_value B on A.seqfeature_id=B.seqfeature_id ' \
           'inner join term as C on B.term_id=C.term_id where C.name="product"' % (biodb, filter)
    #print sql3
    seqfeature_id2product = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))

    sql4 = 'select A.seqfeature_id,start_pos,end_pos,strand ' \
           ' from (select t2.bioentry_id,taxon_id,accession,seqfeature_id,type_term_id,t4.name ' \
           ' from biodatabase as t1 ' \
           ' inner join bioentry as t2 on t1.biodatabase_id =t2.biodatabase_id ' \
           ' inner join seqfeature as t3 on t2.bioentry_id=t3.bioentry_id ' \
           ' inner join term t4 on t3.type_term_id=t4.term_id where t1.name="%s" and t4.name in (%s))A ' \
           ' inner join location B on A.seqfeature_id=B.seqfeature_id' % (biodb, filter)

    seqfeature_id2location = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4,))
    #print sql4

    for seqfeature_data in non_protein_coding_data:

        bioentry_id, taxon_id, accession, seqfeature_id, type_term_id, name = seqfeature_data
        #print type(seqfeature_id)
        try:
            locus_tag = seqfeature_id2locus_tag[str(seqfeature_id)]
        except:
            locus_tag = '-'
        try:
            product = seqfeature_id2product[str(seqfeature_id)]
        except:
            product = '-'
        start, stop, strand = seqfeature_id2location[seqfeature_id]

        sql = 'insert into non_protein_coding_locus values (%s,%s,%s,"%s","%s","%s","%s",%s,%s,%s)' % (seqfeature_id,
                                                                                                       taxon_id,
                                                                                                       bioentry_id,
                                                                                                       accession,
                                                                                                       name,
                                                                                                       locus_tag,
                                                                                                       product,
                                                                                                       start,
                                                                                                       stop,
                                                                                                       strand)

        server.adaptor.execute(sql,)
        server.adaptor.commit()


def get_accession_list_from_taxon_id(server, biodatabase_name, taxon_id):
    sql = "select accession from bioentry" \
          " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s" \
          " where taxon_id = %s"
    result = server.adaptor.execute_and_fetchall(sql, (biodatabase_name, taxon_id))
    return [i[0] for i in result]

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import biosql_own_sql_tables
    from chlamdb.biosqldb import get_locus2seqfeature_table
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", '--mcl', type=str, help="mcl file")
    parser.add_argument("-d", '--db_name', type=str, help="db name")

    args = parser.parse_args()

    server, db = manipulate_biosqldb.load_db(args.db_name)
    asset_path = "/home/trestan/work/dev/django/chlamydia/assets"


    print("parsing orthofinder file")
    locus_tag2orthogroup_id, \
    orthomcl_groups2locus_tag_list, \
    genome_orthomcl_code2proteins, \
    protein_id2genome_ortho_mcl_code = parse_orthomcl_output(args.mcl,
                                                             True)

    print("get locus_tag2seqfeature_id")
    locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, args.db_name)
    locus_tag2seqfeature_id_CDS = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, args.db_name, all=False)

    print("number of groups:", len(orthomcl_groups2locus_tag_list))
    print("number of locus in locus_tag2orthogroup_id:", len(locus_tag2orthogroup_id))
    print("number of locus in locus_tag2seqfeature_id:", len(locus_tag2seqfeature_id))
    print("number of locus in locus_tag2seqfeature_id_CDS:", len(locus_tag2seqfeature_id_CDS))

    
    print("adding orthogroup to seqfeature_qualifier_values")   
    add_orthogroup_to_seq(server,
                            locus_tag2orthogroup_id,
                            locus_tag2seqfeature_id)
    
    
    print("Get orthology matrix merging plasmid")
    orthogroup2detailed_count = get_orthology_matrix_merging_plasmids_biosqldb(server, args.db_name)

    
    print("creating orthology table")
    print("number of groups", len(orthogroup2detailed_count))
    create_orthology_mysql_table(server, orthogroup2detailed_count, args.db_name)
    
    
    print("get locus_tag2taxon_id dictionnary...")
    locus_tag2genome_taxon_id = manipulate_biosqldb.locus_tag2genome_taxon_id(server, args.db_name)


    
    print('creating locustag2seqfature_id table')
    get_locus2seqfeature_table.create_locus_tag2seqfeature_table(args.db_name,
                                                                 locus_tag2seqfeature_id,
                                                                 locus_tag2genome_taxon_id)
    
    
    print("get protein_id2seqfeature_id")
    protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, args.db_name)

    print("get locus_tag2accession dictionnary...")
    locus_tag2accession = manipulate_biosqldb.locus_tag2accession(server, args.db_name)

    print("get protein_id2accession dictionnary...")
    protein_id2accession = manipulate_biosqldb.protein_id2accession(server, args.db_name)

    print("getting location")
    seqfeature_id2seqfeature_location = manipulate_biosqldb.seqfeature_id2feature_location_dico(server, args.db_name)

    print("getting seqfeature_id2protein_id")
    seqfeature_id2protein_id = manipulate_biosqldb.seqfeature_id2protein_id_dico(server, args.db_name)

    print("getting seqfeature_id2bioentry_id")
    seqfeature_id2bioentry_id = manipulate_biosqldb.seqfeature_id2bioentry_id_dico(server, args.db_name)

    print("getting seqfeature_id2gene")
    seqfeature_id2gene = manipulate_biosqldb.seqfeature_id2gene_dico(server, args.db_name)

    print("getting seqfeature_id2product")
    seqfeature_id2product = manipulate_biosqldb.seqfeature_id2product_dico(server, args.db_name)

    print ("getting seqfeature_id2feature_type_id")
    seqfeature_id2feature_type_id = manipulate_biosqldb.seqfeature_id2feature_type_id(server, args.db_name)

    print ("getting seqfeature_id2translation")
    seqfeature_id2translation = manipulate_biosqldb.seqfeature_id2translation_dico(server, args.db_name)

    print ("getting seqfeature_id2organism")
    seqfeature_id2organism = manipulate_biosqldb.seqfeature_id2organism_dico(server, args.db_name)

    print("get pseudogene feature list")
    pseudogene_feature_list = manipulate_biosqldb.pseudogene_feature_list(server, args.db_name)

    print('getting group2group size')
    group2group_size = get_all_orthogroup_size(server, args.db_name)

    print("getting seqfeature_id2locus_tag")
    seqfeature_id2locus_tag = manipulate_biosqldb.seqfeature_id2locus_tag_dico(server, args.db_name)

    print('getting group2family size')
    group2family_size = get_family_size(server, args.db_name)

    create_orthology_tables(server)
    
    print("creating annotation tables (seqfeature_id2locus, seqfeature_id2cds, seqfeature_id2rrna)")
    create_annotation_seqfeature_id2locus_and_cds(server,
                                                  locus_tag2seqfeature_id,
                                                  locus_tag2genome_taxon_id,
                                                  seqfeature_id2protein_id,
                                                  seqfeature_id2gene,
                                                  seqfeature_id2product,
                                                  seqfeature_id2translation,
                                                  seqfeature_id2bioentry_id,
                                                  seqfeature_id2seqfeature_location,
                                                  seqfeature_id2feature_type_id,
                                                  pseudogene_feature_list)
    
    
    print("creating orthology_orthogroup table 1")
    create_orthology_orthogroup(server,
                                group2group_size,
                                group2family_size,
                                orthomcl_groups2locus_tag_list)
    
    print("creating eqfeature_id2orthogroup table 1")
    orthology_seqfeature_id2orthogroup(server,
                                       orthomcl_groups2locus_tag_list,
                                       locus_tag2seqfeature_id)
    

    print ("creating protein_id2taxon_id dictionnary...")
    protein_id2genome_taxon_id = manipulate_biosqldb.protein_id2genome_taxon_id(server, args.db_name)

    print ("creating legacy orthology detail table...")
    create_orthogroup_table_legacy(server,
                                   args.db_name,
                                   orthomcl_groups2locus_tag_list,
                                   locus_tag2seqfeature_id,
                                   protein_id2seqfeature_id,
                                   locus_tag2genome_taxon_id,
                                   protein_id2genome_taxon_id,
                                   locus_tag2accession,
                                   protein_id2accession,
                                   seqfeature_id2locus_tag,
                                   seqfeature_id2protein_id,
                                   seqfeature_id2gene,
                                   seqfeature_id2product,
                                   seqfeature_id2translation,
                                   seqfeature_id2organism,
                                   seqfeature_id2seqfeature_location,
                                   group2group_size,
                                   group2family_size)
    
    # add indexes
    print ("Indexing columns...")
    create_indexes(server)

    # update config
    manipulate_biosqldb.update_config_table(args.db_name, "orthology_data")