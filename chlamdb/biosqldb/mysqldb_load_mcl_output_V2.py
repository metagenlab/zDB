#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-



from parse_mclOtput import parse_orthomcl_output
from chlamdb.biosqldb import manipulate_biosqldb
import os
import shell_command
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
        print "adding orthogroup seqfeature term to", server
        # add orthogroup, ontology term: 2
        sql2 = 'INSERT INTO term (name, ontology_id) values ("orthogroup", 2);'
        server.adaptor.execute(sql2)
        server.adaptor.commit()
    return server.adaptor.execute_and_fetchall(sql1)[0][0]

def add_orthogroup_to_seq(server, protein_id2orthogroup, protein_id2seqfeature_id_dico, locus_tag2seqfeature_id_dico):
    #| seqfeature_id | term_id | rank | value   |
    rank = 1
    term_id = add_orthogroup_term(server)
    for protein_id in protein_id2orthogroup.keys():
        try:
            seqfeature_id = protein_id2seqfeature_id_dico[protein_id]
        except:
            seqfeature_id = locus_tag2seqfeature_id_dico[protein_id]
        group = protein_id2orthogroup[protein_id]
        sql = 'INSERT INTO seqfeature_qualifier_value (seqfeature_id, term_id, rank, value) values (%s, %s, %s, "%s");' % (seqfeature_id,
                                                                                                                           term_id,
                                                                                                                           rank,
                                                                                                                           group)
        try:
            server.adaptor.execute(sql)
        except:
            print 'group %s already inserted?' % group
            print 'updating...'
            #try:
            sql2 = 'UPDATE seqfeature_qualifier_value SET seqfeature_qualifier_value.value="%s" where seqfeature_id=%s and term_id=%s' % (group,
                                                                                                                seqfeature_id,
                                                                                                                term_id)
            #print sql2
            server.adaptor.execute(sql2)
            #except:
            #    print 'could not insert group %s' % group


            print sql

    server.adaptor.commit()
        #select term_id from term where name = "orthogroup";







def create_orthogroup_table(server, biodatabase_name,
                            orthomcl_groups2locus_tag_list):
    import re

    sql = 'CREATE TABLE if not EXISTS orthologous_groups_%s(orthogroup VARCHAR(100) NOT NULL, ' \
          ' protein_id INT, ' \
          ' orthogroup_size INT,' \
          ' n_genomes INT, INDEX protein_id(protein_id), INDEX orthogroup(orthogroup))' % biodatabase_name

    server.adaptor.execute(sql)


    for i in range(0, len(orthomcl_groups2locus_tag_list.keys())):
        group = "group_%s" % i

        locus_tag_list = orthomcl_groups2locus_tag_list[group]
        orthogroup_size = len(locus_tag_list)

        filter = '"' + '","'.join(locus_tag_list) + '"'
        sql = 'select * from feature_tables.genomes_cds_%s as t1 inner join feature_tables.cds_accessions_%s as t2 on ' \
              ' t1.prot_primary_id=t2.prot_primary_id where id_type="locus_tag" and t2.id_name in (%s) ' \
              ' group by taxon_id' % (biodatabase_name,biodatabase_name,
                                      filter)
        print sql
        n_genomes = len(server.adaptor.execute_and_fetchall(sql,))

        # iterate all proteins from orthomcl file
        # add data to the summary table othology_detail...
        print group, 'n genomes', n_genomes
        for protein_id in locus_tag_list:

            sql = 'INSERT INTO orthologous_groups_%s(orthogroup, protein_id,orthogroup_size, n_genomes) ' \
                  ' values ("%s", "%s", %s, %s);' % (biodatabase_name,
                                                       group,
                                                       protein_id,
                                                       orthogroup_size,
                                                       n_genomes)
            try:
                server.adaptor.execute(sql)
                server.adaptor.commit()
            except:
                print 'problem with:'
                print sql


def get_all_orthogroup_size(server, biodatabase_name):
    """
    return a dictonary with orthogroup size"
    """

    sql = ' select seqfeature_qualifier_value.value, COUNT(*) from seqfeature_qualifier_value' \
          ' inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = "orthogroup"' \
          ' inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id' \
          ' inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id' \
          ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "%s"' \
          ' group by seqfeature_qualifier_value.value' % biodatabase_name

    print sql

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
            print sql, key
            description = server.adaptor.execute_and_fetchall(sql, key)[0][0]
        
        print all_taxons
        all_data = {}
        for taxon in all_taxons:
            sql = 'select orthogroup,`%s` from orthology_%s where `%s` >0' % (taxon, biodatabase_name, taxon)
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
        print n
        if n%2 ==0:
            n+=1
        fig, axes = plt.subplots(nrows=n, ncols=2)
        for key, location in zip(complete_data.keys(), list(itertools.product(range(n),repeat=2))):
            serie = complete_data[key]
            sql = 'select description from bioentry where taxon_id=%s and description not like "%%%%plasmid%%%%" limit 1;'
            print sql, key
            description = server.adaptor.execute_and_fetchall(sql, key)[0][0]
            print description

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
        sql='select orthogroup, count(*) from orthology_detail group by orthogroup' % biodatabase_name

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
        print "ID!!!"
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
        print genome
        dico = manipulate_biosqldb.bioentry_name2orthoup_size(server, biodatabase_name, genome)
        print dico



def get_orthology_matrix_merging_plasmids_biosqldb(server, biodatabase_name):

    all_orthogroups = get_all_orthogroup_size(server, biodatabase_name)
    #server, db = manipulate_biosqldb.load_db(biodatabase_name)
    #sql='select orthogroup, count(*) from orthology_detail group by orthogroup' % biodatabase_name

    # based on the new orthology detail table
    #all_orthogroups = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


    n_groups = len(all_orthogroups.keys())
    #print "ngroups", n_groups


    all_taxons = manipulate_biosqldb.get_taxon_id_list(server, biodatabase_name)
    #print all_taxons
    print 'get dico ortho count'
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
    print 'count ok'
    #print detailed_orthology_count
    return detailed_orthology_count



def get_orthology_matrix_merging_plasmids_own_tables(server, biodatabase_name):

    #all_orthogroups = get_all_orthogroup_size(server, biodatabase_name)
    server, db = manipulate_biosqldb.load_db(biodatabase_name)
    sql='select orthogroup, count(*) from orthology_detail group by orthogroup' % biodatabase_name

    # based on the new orthology detail table
    all_orthogroups = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


    n_groups = len(all_orthogroups.keys())
    #print "ngroups", n_groups

    
    all_taxons = manipulate_biosqldb.get_taxon_id_list(server, biodatabase_name)
    #print all_taxons
    print 'get dico ortho count'
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
    print 'count ok'
    #print detailed_orthology_count
    return detailed_orthology_count



def create_orthology_mysql_table(server, orthogroup2detailed_count, biodatabase_name):

    """
    CREATE TABLE table_name (column_name column_type);    
    """
    print 
    taxon_id_list = orthogroup2detailed_count[orthogroup2detailed_count.keys()[0]].keys()
    
    sql = "CREATE TABLE comparative_tables_orthology(orthogroup VARCHAR(100) NOT NULL" % biodatabase_name
    for i in taxon_id_list:
        sql+=" ,`%s` INT" % i
    sql+=")"
    server.adaptor.execute(sql)

    for group in orthogroup2detailed_count.keys():
        values='"%s", ' % group
        columns='orthogroup, '
        taxons = orthogroup2detailed_count[group].keys()
        for i in range(0, len(taxons)-1):
            values += " %s," % orthogroup2detailed_count[group][taxons[i]]
            columns+="`%s`, " % taxons[i]
        values += "%s" % orthogroup2detailed_count[group][taxons[-1]]
        columns += "`%s`" % taxons[-1]
        sql = "INSERT INTO comparative_tables_orthology (%s) values (%s);" %  (biodatabase_name, columns, values)
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





    
def get_one_group_data(group_list, biodatabase_name, out_dir):#, out_q):
    server = manipulate_biosqldb.load_db()
    for group in group_list:
        print "group", group, "biodb", biodatabase_name
        seqfeature_ids = manipulate_biosqldb.orthogroup_id2seqfeature_id_list(server, group, biodatabase_name)
        # do not consider singletons
        one_fasta = ""
        if len(seqfeature_ids) > 0:
            for seqfeature_id in seqfeature_ids:
                # get data
                values = manipulate_biosqldb.seqfeature_id2seqfeature_qualifier_values(server, seqfeature_id, biodatabase_name)
                try:
                    one_fasta += ">%s | %s, %s | %s\n%s\n" % (values['locus_tag'],  values['gene'], values['product'], values['protein_id'], values['translation'])
                except:
                    try:
                        one_fasta += ">%s | %s, %s\n%s\n" % (values['locus_tag'], values['gene'], values['product'], values['translation'])
                    except:
                        try:
                            one_fasta += ">%s | %s\n%s\n" % (values['locus_tag'], values['product'], values['translation'])
                        except:
                            try:
                                #print values
                                one_fasta += ">%s \n%s\n" % (values['locus_tag'], values['translation'])
                            except:
                                one_fasta += ">%s \n%s\n" % (values['protein_id'], values['translation'])
            out_name = os.path.join(out_dir, "%s.txt" % group)
            f = open(out_name, "w")
            f.write(one_fasta)
        f.close()
        #out_q.put(group)

    

def get_all_orthogroup_protein_fasta(server, biodatabase_name, out_dir):

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodatabase_name)
    sql = 'select orthogroup from orthology_detail group by orthogroup' % biodatabase_name


    all_groups = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)] # get_all_orthogroup_size(server, biodatabase_name).keys()
    #all_groups2 =  list(get_all_orthogroup_size(server, biodatabase_name).keys())
    #print all_groups2[0:5], all_groups[0:5]
    #print 'diff', list(set(all_groups2) - set(all_groups))

    #out_q = Queue()
    n_cpu = 8
    n_poc_per_list = int(numpy.ceil(len(all_groups)/float(n_cpu)))
    query_lists = _chunks(all_groups, n_poc_per_list)
    procs = []
    for one_list in query_lists:
        proc = Process(target=get_one_group_data, args=(one_list, biodatabase_name, out_dir))#, out_q))
        procs.append(proc)
        proc.start()

        
    print "join proc"
    time.sleep(5)
    for proc in procs:
        proc.join()
        

def get_all_orthogroup_protein_fasta2(server, biodatabase_name, out_dir):

    all_groups = get_all_orthogroup_size(server, biodatabase_name).keys()
    
    for group in all_groups:
        print "group", group
        seqfeature_ids = manipulate_biosqldb.orthogroup_id2seqfeature_id_list(server, group, biodatabase_name)
        # do not consider singletons
        one_fasta = ""
        if len(seqfeature_ids) > 0:
            for seqfeature_id in seqfeature_ids:
                # get data
                values = manipulate_biosqldb.seqfeature_id2seqfeature_qualifier_values(server, seqfeature_id, biodatabase_name)
                try:
                    one_fasta += ">%s | %s, %s | %s\n%s\n" % (values['locus_tag'],  values['gene'], values['product'], values['protein_id'], values['translation'])
                except:
                    try:
                        one_fasta += ">%s | %s, %s\n%s\n" % (values['locus_tag'], values['gene'], values['product'], values['translation'])
                    except:
                        try:
                            one_fasta += ">%s | %s\n%s\n" % (values['locus_tag'], values['product'], values['translation'])
                        except:
                            try:
                                #print values
                                one_fasta += ">%s \n%s\n" % (values['locus_tag'], values['translation'])
                            except:
                                one_fasta += ">%s \n%s\n" % (values['protein_id'], values['translation'])
            out_name = os.path.join(out_dir, "%s.txt" % group)
            f = open(out_name, "w")
            f.write(one_fasta)
        f.close() 









def get_one_group_data_taxon(group_list, locus_or_protein_id2taxon_id, biodatabase_name, out_dir):#, out_q):
    server = manipulate_biosqldb.load_db()
    for group in group_list:
        print "group", group
        seqfeature_ids = manipulate_biosqldb.orthogroup_id2seqfeature_id_list(server, group, biodatabase_name)
        # do not consider singletons
        one_fasta = ""
        if len(seqfeature_ids) > 0:
            for seqfeature_id in seqfeature_ids:
                # get data
                values = manipulate_biosqldb.seqfeature_id2seqfeature_qualifier_values(server, seqfeature_id, biodatabase_name)
                #print values
                try:
                    taxon = locus_or_protein_id2taxon_id[values['locus_tag']]
                    one_fasta += ">%s \n%s\n" % (taxon, values['translation'])
                except:
                    taxon = locus_or_protein_id2taxon_id[values['protein_id']]
                    one_fasta += ">%s \n%s\n" % (taxon, values['translation'])
            out_name = os.path.join(out_dir, "%s.txt" % group)
            f = open(out_name, "w")
            f.write(one_fasta)
        #out_q.put(group)


        
def get_all_orthogroup_protein_fasta_by_taxon(server, biodatabase_name, out_dir):

    all_groups = get_all_orthogroup_size(server, biodatabase_name).keys()
    locus_or_protein_id2taxon_id = manipulate_biosqldb.locus_or_protein_id2taxon_id(server, biodatabase_name)

    #out_q = Queue()
    n_cpu = 8
    n_poc_per_list = int(numpy.ceil(len(all_groups)/float(n_cpu)))
    query_lists = _chunks(all_groups, n_poc_per_list)
    procs = []
    
    for one_list in query_lists:
        proc = Process(target=get_one_group_data_taxon, args=(one_list, locus_or_protein_id2taxon_id, biodatabase_name, out_dir))#, out_q))
        procs.append(proc)
        proc.start()
        
    print "join proc"
    time.sleep(5)
    for proc in procs:
        proc.join()
    

'''            
def circos(server, biodatabase_name, reference_taxon_id):

    
    
    import gbk2circos

    accessions = get_accession_list_from_taxon_id(server, biodatabase_name, reference_taxon_id)
    #for accession in accessions:
    print accessions
    # ajouter plasmide SneZ
    
    record_list = [ gbk2circos.Record(db.lookup(accession=accession)) for accession in accessions]
    #print record_list
    circos_files, taxon_id2description = gbk2circos.orthology_circos_files(server, record_list, reference_taxon_id, biodatabase_name)
    print circos_files
    print "taxon_id2description", taxon_id2description
    corcos = gbk2circos.Circos_config(circos_files["contigs"])

    # add GC group_size
    
    corcos.add_plot(circos_files["GC"], type="line", r0="1.01r", r1="1.1r", color="green", fill_color="vlgreen", thickness = "2p", z = 1, rules ="")
    corcos.add_plot(circos_files["orthogroups"], type="line", r0="1.12r", r1= "1.22r", color="black", fill_color="red", thickness = "2p", z = 1, rules ="")
        

    # add plus minus genes    
    corcos.add_highlight(circos_files["plus"], fill_color="violet", r1="0.98r", r0="0.95r")
    corcos.add_highlight(circos_files["minus"], fill_color="orange", r1="0.95r", r0="0.92r")

    # add presence/absence of orthologs

    r1 = 0.90
    r0 = 0.87
    for orthofile in circos_files["genomes"]:
     print orthofile
     corcos.add_highlight(orthofile, fill_color="blue", r1= "%sr" % r1, r0= "%sr" % r0)
     r1 = r1-0.05
     r0 = r0-0.05

    accessions_name = ""
    for accession in accessions:
        accessions_name+="_" + accession

    config_file = "circos_config%s.txt" % accessions_name
    t = open(config_file, "w")
    t.write(corcos.get_file())
    
    print taxon_id2description
    return (config_file, "circos" + accessions_name)



    
def circos_draft(server, biodatabase_name, reference_taxon_id, fasta_draft_reference):

    
    
    import gbk2circos


    draft_contigs = circos_fasta_draft(fasta_draft_reference)
    print "draft_contigs", draft_contigs
    accessions = get_accession_list_from_taxon_id(server, biodatabase_name, reference_taxon_id)
    #for accession in accessions:
    print accessions
    # ajouter plasmide SneZ
    
    record_list = [ gbk2circos.Record(db.lookup(accession=accession)) for accession in accessions]
    #print record_list
    circos_files, taxon_id2description = gbk2circos.orthology_circos_files(server, record_list, reference_taxon_id, biodatabase_name, draft_data = draft_contigs)
    print circos_files
    corcos = gbk2circos.Circos_config(circos_files["contigs"])

    # add GC group_size
    
    corcos.add_plot(circos_files["GC"], type="line", r0="1.01r", r1="1.1r", color="green", fill_color="vlgreen", thickness = "2p", z = 1, rules ="")
    corcos.add_plot(circos_files["orthogroups"], type="line", r0="1.12r", r1= "1.22r", color="black", fill_color="red", thickness = "2p", z = 1, rules ="")
        

    # add plus minus genes    
    corcos.add_highlight(circos_files["plus"], fill_color="grey_a1", r1="0.98r", r0="0.95r")
    corcos.add_highlight(circos_files["minus"], fill_color="grey_a1", r1="0.95r", r0="0.92r")

    # add presence/absence of orthologs

    r1 = 0.90
    r0 = 0.89
    for orthofile in circos_files["genomes"]:
     print orthofile
     corcos.add_highlight(orthofile, fill_color="blue", r1= "%sr" % r1, r0= "%sr" % r0)
     r1 = r1-0.02
     r0 = r0-0.02

    accessions_name = ""
    for accession in accessions:
        accessions_name+="_" + accession

    config_file = "circos_config%s.txt" % accessions_name
    t = open(config_file, "w")
    t.write(corcos.get_file())
    
    print taxon_id2description
    return (config_file, "circos" + accessions_name)
'''

    
def get_conserved_core_groups(server, biodatabase_name):
    all_taxons_id = manipulate_biosqldb.get_taxon_id_list(server, biodatabase_name)


    sql_taxons = ""
    for i in range(0, len(all_taxons_id)-1):
        sql_taxons += ' `%s` = 1 and' % all_taxons_id[i]
    sql_taxons += ' `%s` = 1' % all_taxons_id[-1]
        
    sql = "select orthogroup from orthology_%s where %s;" % (biodatabase_name, sql_taxons)
    #print sql
    result = server.adaptor.execute_and_fetchall(sql,)
    
    return [i[0] for i in result]


def get_cds_sequences(record, locus_tag2taxon_id):
    from Bio.SeqRecord import SeqRecord
    cds = {}
    for i in record.features:
        if i.type == "CDS":
            try:
                taxon_id = locus_tag2taxon_id[str(i.qualifiers['locus_tag'][0])]
                #print str(i.qualifiers['locus_tag'][0])
                if str(i.qualifiers['locus_tag'][0]) == "HMPREF0772_11479":
                    print i
                    print len(i)
                    print i.location
                    print i.extract(record.seq)
                    print len(i.extract(record.seq))
                    #import sys
                    #sys.exit()
                #if i.location.strand == -1:
                #    print "reverse!"
                #seq = i.extract(record.seq)[::-1]
                #else:
                #    print "not reverse!"
                seq = i.extract(record.seq)
                cds[str(i.qualifiers['locus_tag'][0])] = SeqRecord(seq, id=str(taxon_id), description=i.qualifiers['locus_tag'][0], name=i.qualifiers['locus_tag'][0])
            except:
                #sys.exit()
                print i
    return cds


def get_group_nucl_fasta(biodb_name, locus_tag2nucl_sequence, path, group_list):
    from Bio import SeqIO

    for group in group_list:
        server = manipulate_biosqldb.load_db()
        locus_list = manipulate_biosqldb.orthogroup_id2locus_tag_list(server, group, biodb_name)
        seqs = []

        for locus in locus_list:
            #print locus
            seq = locus_tag2nucl_sequence[str(locus[2])]
            #print seq
            seqs.append(locus_tag2nucl_sequence[str(locus[2])])
        with open(os.path.join(path, group + "_nucl.txt"), "w") as f:
            SeqIO.write(seqs, f, "fasta")



def locus_tag2nucl_sequence_dict(server, db, biodb_name):
    sql = 'select accession from bioentry' \
          ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "%s"' % biodb_name
    print sql

    accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    locus_or_protein_id2taxon_id = manipulate_biosqldb.locus_or_protein_id2taxon_id(server, biodb_name)

    locus_tag2nucl_sequence = {}
    for accession in accessions:
        print accession
        record_data = db.lookup(accession=accession)
        sequences = get_cds_sequences(record_data, locus_or_protein_id2taxon_id)
        locus_tag2nucl_sequence.update(sequences)
    return locus_tag2nucl_sequence

def locus_tag2aa_sequence_dict(server, db, biodb_name):

    sql = 'select locus_tag, translation from orthology_detail' % biodb_name

    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

def orthogroup2nucleotide_seq_list(locus_tag2nucl_sequence, group, biodb_name):
        server = manipulate_biosqldb.load_db()
        locus_list = manipulate_biosqldb.orthogroup_id2locus_tag_list(server, group, biodb_name)
        seqs = []
        for locus in locus_list:
            seq = locus_tag2nucl_sequence[str(locus[2])]
            seqs.append(locus_tag2nucl_sequence[str(locus[2])])
        return seqs


def get_nucleotide_core_fasta(server, db, biodb_name, path):
    from Bio import SeqIO

    all_groups = get_all_orthogroup_size(server, biodb_name).keys()
    #core_ortho = get_conserved_core_groups(server, biodb_name)

    my_locus_tag2nucl_sequence = locus_tag2nucl_sequence_dict(server, db, biodb_name)

    for group in all_groups:
        seqs = orthogroup2nucleotide_seq_list(my_locus_tag2nucl_sequence, group, biodb_name)
        with open(os.path.join(path,  group + "_nucl.txt"), "w") as f:
            SeqIO.write(seqs, f, "fasta")




    '''
    n_cpu = 8
    n_poc_per_list = int(numpy.ceil(len(core_ortho)/float(n_cpu)))
    query_lists = _chunks(core_ortho, n_poc_per_list)
    procs = []

    for one_list in query_lists:
        proc = Process(target=get_group_nucl_fasta, args=(biodb_name, locus_tag2nucl_sequence, path, one_list))#, out_q))
        procs.append(proc)
        proc.start()

    print "join proc"
    time.sleep(5)
    for proc in procs:
        proc.join()

    '''


    
def get_accession_list_from_taxon_id(server, biodatabase_name, taxon_id):
    sql = "select accession from bioentry" \
          " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s" \
          " where taxon_id = %s" 
    result = server.adaptor.execute_and_fetchall(sql, (biodatabase_name, taxon_id))
    return [i[0] for i in result]

if __name__ == '__main__':
    import argparse
    import biosql_own_sql_tables

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", '--mcl',type=str,help="mcl file")
    parser.add_argument("-d", '--db_name', type=str,help="db name")
    parser.add_argument("-f", '--fasta_draft', type=str,help="draft reference genome")
    parser.add_argument("-p", '--asset_path', type=str,help="asset path")
    parser.add_argument("-t", '--locus_tag2protein_id_files', type=str,help="locus_tag2protein_id files", nargs='+')
    parser.add_argument("-s", '--get_sequences', action="store_true", help="Create aa and nucleotide alignment directories in asset path")
    parser.add_argument("-c", '--core_groups_path', type=str, help="Path to save core orthogroup fasta. Taxon id as header for concatenation.", default=None)
    parser.add_argument("-o", '--orthofinder', action="store_true",help="orthofinder input file (and not orthomcl)")
    args = parser.parse_args()
    
    server, db = manipulate_biosqldb.load_db(args.db_name)
    asset_path = "/home/trestan/work/dev/django/chlamydia/assets"


    if not args.get_sequences and not args.core_groups_path:

        print "creating locus_tag2seqfeature_id"
        locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, args.db_name)


        print "creating protein_id2seqfeature_id"
        protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, args.db_name)

        print "creating locus_tag2taxon_id dictionnary..."
        locus_tag2genome_taxon_id = manipulate_biosqldb.locus_tag2genome_taxon_id(server, args.db_name)
        print "creating protein_id2taxon_id dictionnary..."
        protein_id2genome_taxon_id = manipulate_biosqldb.protein_id2genome_taxon_id(server, args.db_name)

        print "creating locus_tag2accession dictionnary..."
        locus_tag2accession = manipulate_biosqldb.locus_tag2accession(server, args.db_name)
        print "creating protein_id2accession dictionnary..."
        protein_id2accession = manipulate_biosqldb.protein_id2accession(server, args.db_name)

        print "getting location"
        seqfeature_id2seqfeature_location = manipulate_biosqldb.seqfeature_id2feature_location_dico(server, args.db_name)

        print "getting seqfeature_id2locus_tag"
        seqfeature_id2locus_tag = manipulate_biosqldb.seqfeature_id2locus_tag_dico(server, args.db_name)

        print "getting seqfeature_id2protein_id"
        seqfeature_id2protein_id = manipulate_biosqldb.seqfeature_id2protein_id_dico(server, args.db_name)

        print "parsing mcl file"
        protein_id2orthogroup_id, \
        orthomcl_groups2locus_tag_list, \
        genome_orthomcl_code2proteins, \
        protein_id2genome_ortho_mcl_code = parse_orthomcl_output(args.mcl,
                                                                 args.orthofinder)

        print "getting seqfeature_id2gene"
        seqfeature_id2gene = manipulate_biosqldb.seqfeature_id2gene_dico(server, args.db_name)

        print "getting seqfeature_id2product"
        seqfeature_id2product = manipulate_biosqldb.seqfeature_id2product_dico(server, args.db_name)

        print "getting seqfeature_id2translation"
        seqfeature_id2translation = manipulate_biosqldb.seqfeature_id2translation_dico(server, args.db_name)

        print "getting seqfeature_id2organism"

        seqfeature_id2organism = manipulate_biosqldb.seqfeature_id2organism_dico(server, args.db_name)

        '''
        print "adding orthogroup to seqfeature_qualifier_values"
        add_orthogroup_to_seq(server, protein_id2orthogroup_id, protein_id2seqfeature_id, locus_tag2seqfeature_id)
        '''

        print "creating orthology table merging plasmid"
        orthogroup2detailed_count = get_orthology_matrix_merging_plasmids_biosqldb(server, args.db_name)
        print "orthogroup2detailed_count", orthogroup2detailed_count

        print "creating orthology table"
        create_orthology_mysql_table(server, orthogroup2detailed_count, args.db_name)


        #print "creating orthology table 1"
        #create_orthogroup_table(server, args.db_name,
        #                        orthomcl_groups2locus_tag_list)


        print 'adding orthogroup to interpro table'
        biosql_own_sql_tables.add_orthogroup_to_interpro_table(args.db_name)
        
        print "plotting orthogroup_size"
        plot_orthogroup_size_distrib(server, args.db_name)

        #print "writing fasta files"
        #ortho_table = "orthology_%s" % args.db_name
        #all_taxon_ids = manipulate_biosqldb.get_column_names(server, ortho_table)[1:]
        
    if args.get_sequences and not args.core_groups_path:

        if not os.path.exists(os.path.join(asset_path, "%s_fasta/" % args.db_name)):
            os.makedirs(os.path.join(asset_path, "%s_fasta/" % args.db_name))

        get_all_orthogroup_protein_fasta(server, args.db_name, os.path.join(asset_path, "%s_fasta/" % args.db_name))

        if not os.path.exists(os.path.join(asset_path, "%s_fasta_by_taxons/" % args.db_name)):
            os.makedirs(os.path.join(asset_path, "%s_fasta_by_taxons/" % args.db_name))

        get_all_orthogroup_protein_fasta_by_taxon(server, args.db_name, os.path.join(asset_path, "%s_fasta_by_taxons/" % args.db_name))

        if not os.path.exists(os.path.join(asset_path, "%s_fasta_core/" % args.db_name)):
            os.makedirs(os.path.join(asset_path, "%s_fasta_core/" % args.db_name))
        
        print 'getting single copy core groups'    
        core_ortho = get_conserved_core_groups(server, args.db_name)
        import shutil
        for group in core_ortho:
            shutil.copy(os.path.join(asset_path, "%s_fasta_by_taxons/%s.txt" % (args.db_name, group)),
                    os.path.join(asset_path, "%s_fasta_core/%s.txt" % (args.db_name, group)))

        get_nucleotide_core_fasta(server, db, args.db_name, "./nucl")




    if args.core_groups_path:
        taxon_ids = manipulate_biosqldb.get_taxon_id_list(server, args.db_name)

        sql_include = ''
        if len(taxon_ids) > 0:
            for i in range(0, len(taxon_ids)-1):
                sql_include += ' `%s` = 1 and ' % taxon_ids[i]
            sql_include+='`%s` = 1' % taxon_ids[-1]
        sql ='select orthogroup from orthology_%s where %s' % (args.db_name, sql_include)
        print sql
        group_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
        locus_or_protein_id2taxon_id = manipulate_biosqldb.locus_or_protein_id2taxon_id(server, args.db_name)
        print 'Number of core groups: %s' % len(group_list)
        get_one_group_data_taxon(group_list, locus_or_protein_id2taxon_id, args.db_name, args.core_groups_path)
