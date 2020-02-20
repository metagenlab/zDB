#! /usr/bin/env python


import sys
from Bio import GenBank
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from chlamdb.biosqldb.manipulate_biosqldb import load_db
from chlamdb.biosqldb.manipulate_biosqldb import query_yes_no


def chunks(l, n):
    "return sublists of l of minimum length n (work subdivision for the subprocesing module"
    return [l[i:i+n] for i in range(0, len(l), n)]


def sql_euclidian_dist_orthogroups(biodb, one_list, orthogroup2profile):
    from scipy.spatial.distance import euclidean
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    n = len(one_list)
    for i, one_pair in enumerate(one_list):
        #print ("%s/%s" % (i,n))
        dist = euclidean(orthogroup2profile[one_pair[0]], orthogroup2profile[one_pair[1]])
        if dist <= 2.5:
            sql = 'insert into interactions_phylo_profiles_eucl_dist values ("%s", "%s", %s);' % (one_pair[0],
                                                                                                  one_pair[1],
                                                                                                  dist)
            server.adaptor.execute(sql,)
    server.adaptor.commit()


def sql_jaccard_dist_orthogroups(biodb, one_list, orthogroup2profile):
    from scipy.spatial.distance import jaccard
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    n = len(one_list)
    for i, one_pair in enumerate(one_list):
        #print ("%s/%s" % (i,n))
        dist = jaccard(orthogroup2profile[one_pair[0]], orthogroup2profile[one_pair[1]])
        # scale between 0 and 1
        if dist <= 0.5:
            sql = 'insert into interactions_phylo_profiles_jac_dist values ("%s", "%s", %s);' % (one_pair[0],
                                                                                                 one_pair[1],
                                                                                                 dist)
            server.adaptor.execute(sql,)
    server.adaptor.commit()


def merge_dataframe_columns(dataframe, columns_clusters_dict):
    import pandas

    # create empty dataframe of size n cluster vs n orthogroups
    new_data_frame = pandas.DataFrame(index=dataframe.index.values, columns=columns_clusters_dict.keys())

    # sum row for each cluster of columns (taxons)
    for cluster in columns_clusters_dict:
        new_data_frame[cluster] = dataframe[columns_clusters_dict[cluster]].sum(axis=1)

    return new_data_frame


def create_ortogroup_dist_table(biodb, table_name='phylo_profiles_eucl_dist'):
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)

    sql_profiles_table = 'CREATE TABLE IF NOT EXISTS interactions_%s (group_1 varchar(100), ' \
                         ' group_2 varchar(100), euclidian_dist FLOAT, INDEX group_1 (group_1), INDEX group_2 (group_2))' % (table_name)
                                        
    try:
        print (sql_profiles_table)
        server.adaptor.execute(sql_profiles_table)
    except:
        print ('problem creating the sql table')


def get_reduced_orthogroup_matrix(biodb, jaccard=True):
    '''
    # get orthology matrix
    # perform hierarchical clustering of genome based on median protein identity
    # return reduced matrix
    '''

    import pandas
    import numpy
    import os
    import rpy2.robjects.numpy2ri
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri as numpy2ri
    from chlamdb.biosqldb import manipulate_biosqldb

    rpy2.robjects.numpy2ri.activate()

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select * from comparative_tables_orthology'

    sql2 = 'show columns from comparative_tables_orthology'


    # get matrix as numpy arraw: orthogroups as rows, genomes as columns
    # replace values > 1 by one (profile of presence/absence)
    # possibility: merge closely related data (merge columns)
    #orthogroup_profiles = server.adaptor.execute_and_fetchall(sql,)

    data = numpy.array([list(i) for i in server.adaptor.execute_and_fetchall(sql,)])

    all_cols = [i[0] for i in server.adaptor.execute_and_fetchall(sql2,)]


    # set the persentage of taxon that must have homologs/domain/... => default is 1 (100%)

    orthogroups_df = pandas.DataFrame(data, columns=all_cols)
    groups = orthogroups_df["orthogroup"]
    profiles = orthogroups_df[all_cols[1:]].astype(float)

    #print profiles["131"] + profiles["132"]

    # cluster taxons based on identity of core genome aligment
    # default clustering_ hight of 30

    sql = 'select biodatabase_id from biodatabase where name="%s"'

    biodb_id = server.adaptor.execute_and_fetchall(sql,)[0][0]


    sqlpsw = os.environ['SQLPSW']

    # merge closely related genomes: perform hierarchical clustering based on median protein identity
    # we could use a hard coded cutoff (eg: 95% median protein identity?)
    robjects.r("""
    library("RMySQL")
    library(reshape2)

    con <- dbConnect(MySQL(),
             user="root", password="%s",
             dbname="%s", host="localhost")

    rs1 <- dbSendQuery(con, 'select taxon_1,taxon_2, median_identity from comparative_tables_shared_og_av_id union select taxon_2, taxon_1, median_identity from comparative_tables_shared_og_av_id;')
    pairwise_identity<- dbFetch(rs1, n=-1)
    rs2 <- dbSendQuery(con, 'select taxon_id,description from biosqldb.bioentry where biodatabase_id=%s and description not like "%%plasmid%%";')
    taxon2description<- dbFetch(rs2, n=-1)

    pairwise_identity_matrix <- dcast(pairwise_identity, taxon_1~taxon_2)
    
    rownames(pairwise_identity_matrix) <- pairwise_identity_matrix$taxon_1

    pairwise_identity_matrix<-pairwise_identity_matrix[,2:length(pairwise_identity_matrix)]

    pairwise_dist <- as.dist(100-pairwise_identity_matrix)
    hc <- hclust(pairwise_dist)
    clusterCut <- cutree(hc,h=30)
    taxons <- names(clusterCut)

    """ % (sqlpsw, 
           biodb, 
           biodb_id))

    clusters = robjects.r["clusterCut"]
    taxons = robjects.r["taxons"]

    # build dictionnary of clusters of taxons
    cluster2taxons = {}
    print (type(clusters))
    for clust, tax in zip(clusters, taxons):
        if clust not in cluster2taxons:
            cluster2taxons[clust] = [tax]
        else:
            cluster2taxons[clust].append(tax)

    # build a new DataFrame with merged columns
    print("Merging taxons:")
    print("Cluster\ttaxon")
    for cluster in cluster2taxons:
        for taxon in cluster2taxons[cluster]:
            print("%s\t%s" % (cluster, taxon))
    
    merged_profiles = merge_dataframe_columns(profiles, cluster2taxons)
    # replace multicounts by 1
    merged_profiles[merged_profiles >1] = 1
    
    return groups, merged_profiles


def get_compination_list(groups, merged_profiles):
    '''
    Remove profile with count > 0 in a single species 
    Return 
        all combinations to compare
        dictionnary of profile to compare (group2profile)
    '''
    
    combinations = []

    orthogroup2profile = {}
    for i, group in enumerate(groups):
        orthogroup2profile[group] = merged_profiles.loc[i,]

    filtered_groups = []
    total = len(groups)
    for i, orthogroup in enumerate(groups):
        #if no homologs in other genomes, skip
        if sum(merged_profiles.loc[i,:]) == 1:
            continue
        else:
            filtered_groups.append(orthogroup)

    # get a list of all combination of groups
    # allows then to split the job into 8 process
    total = len(groups)
    for i, orthogroup_1 in enumerate(filtered_groups):
        # if conserved everywhere, skip
        if sum(orthogroup2profile[orthogroup_1]) == len(orthogroup2profile[orthogroup_1]):
            continue
        for y, orthogroup_2 in enumerate(filtered_groups[i:]):# in range(i, total):
            # if group1 == group2, skip
            if orthogroup_1 == orthogroup_2:
                continue
            # if conserved everywhere, skip
            if sum(orthogroup2profile[orthogroup_2]) == len(orthogroup2profile[orthogroup_2]):
                continue
            combinations.append((orthogroup_1, orthogroup_2))
            
    return combinations, orthogroup2profile


def euclidian_dist_orthogroups(biodb, merge_taxons=False, n_cpus=8):

    from chlamdb.biosqldb import manipulate_biosqldb
    from multiprocessing import Process
    import math


    create_ortogroup_dist_table(biodb, 'phylo_profiles_eucl_dist')
    
    groups, merged_profiles = get_reduced_orthogroup_matrix(biodb)

    combinations, orthogroup2profile = get_compination_list(groups, merged_profiles)

    print ('n of pairs:', len(combinations))
    n_poc_per_list = math.ceil(len(combinations)/float(n_cpus))

    # split the jobs into lists (1 per cpu)
    query_lists = chunks(range(0, len(combinations)), int(n_poc_per_list))
    # start processes
    procs = []
    for one_list in query_lists:
        comb_list = [combinations[i] for i in one_list]
        proc = Process(target=sql_euclidian_dist_orthogroups, args=(biodb, comb_list, orthogroup2profile))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()
        


def jaccard_dist_orthogroups(biodb, n_cpus=8):
    
    from chlamdb.biosqldb import manipulate_biosqldb
    from multiprocessing import Process
    import math


    create_ortogroup_dist_table(biodb, 'phylo_profiles_jac_dist')
    
    groups, merged_profiles = get_reduced_orthogroup_matrix(biodb)

    combinations, orthogroup2profile = get_compination_list(groups, merged_profiles)

    print ('n of pairs:', len(combinations))
    n_poc_per_list = math.ceil(len(combinations)/float(n_cpus))

    # split the jobs into lists (1 per cpu)
    query_lists = chunks(range(0, len(combinations)), int(n_poc_per_list))
    # start processes
    procs = []
    for one_list in query_lists:
        comb_list = [combinations[i] for i in one_list]
        proc = Process(target=sql_jaccard_dist_orthogroups, args=(biodb, comb_list, orthogroup2profile))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()



def sql_euclidian_dist_cogs(biodb, one_list, orthogroup2profile):
    import numpy as np
    from scipy.spatial.distance import euclidean
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    n = len(one_list)
    for i, one_pair in enumerate(one_list):
        print ("%s/%s" % (i,n))
        list1 = np.asarray(orthogroup2profile[one_pair[0]])
        list2 = np.asarray(orthogroup2profile[one_pair[1]])
        list1[list1 > 1] = 1
        list2[list2 > 1] = 1

        dist = euclidean(list1, list2)
        if dist <= 2.5:
            sql = 'insert into comparative_tables_cogs_profiles_euclidian_distance values ("%s", "%s", %s);' % (one_pair[0],
                                                                                                                one_pair[1],
                                                                                                                dist)
            server.adaptor.execute(sql,)
            server.adaptor.commit()


def euclidian_dist_cogs(biodb):

    from chlamdb.biosqldb import manipulate_biosqldb
    from multiprocessing import Process
    import math

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_profiles_table = 'CREATE TABLE IF NOT EXISTS comparative_tables_cogs_profiles_euclidian_distance (cog_1 varchar(100), ' \
                         ' cog_2 varchar(100), euclidian_dist FLOAT, INDEX cog_1 (cog_1), INDEX cog_2 (cog_2))'
    try:
        print (sql_profiles_table)
        server.adaptor.execute(sql_profiles_table)
    except:
        print ('problem creating the sql table')

    sql = 'select * from comparative_tables_COG_%s'
    cog2profile = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    combinations = []

    # get a list of all combination of groups
    # allows then to split the job into 8 process
    total = len(cog2profile.keys())
    for i, orthogroup_1 in enumerate(cog2profile.keys()):
        #if i == 1000:
        #    break
        print ("%s/%s" % (i, total))
        if sum(cog2profile[orthogroup_1]) == 1:
            print ('skip!')
            continue
        count = 0
        for n in cog2profile[orthogroup_1]:
            if n>0:
                count+=1
        if count>1:
            print ('range', i, total)
            for y, orthogroup_2 in enumerate(cog2profile.keys()[i:]):# in range(i, total):
                #orthogroup_2 = orthogroup2profile.keys()[y]
                if sum(cog2profile[orthogroup_2]) == 1:
                    continue
                count = 0
                for n in cog2profile[orthogroup_2]:
                    if n>0:
                        count+=1
                if count<=1:
                    continue
                else:
                    combinations.append((orthogroup_1, orthogroup_2))
    n_cpu = 8
    print ('n of pairs:', len(combinations))
    #import sys
    #sys.exit()
    n_poc_per_list = math.ceil(len(combinations)/float(n_cpu))
    #print "n_poc_per_list", n_poc_per_list
    #print range(0, len(combinations))
    # split the jobs into 8 lists
    query_lists = chunks(range(0, len(combinations)), int(n_poc_per_list))
    # start the 8 processes
    procs = []
    for one_list in query_lists:
        comb_list = [combinations[i] for i in one_list]
        #print len(comb_list)
        proc = Process(target=sql_euclidian_dist_cogs, args=(biodb, comb_list, cog2profile))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()


def sql_euclidian_dist_interpro(biodb, one_list, interpro2profile):
    import numpy as np
    from scipy.spatial.distance import euclidean
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    n = len(one_list)
    for i, one_pair in enumerate(one_list):
        #print ("%s/%s" % (i,n))
        list1 = np.asarray(interpro2profile[one_pair[0]])
        list2 = np.asarray(interpro2profile[one_pair[1]])
        list1[list1 > 1] = 1
        list2[list2 > 1] = 1

        dist = euclidean(list1, list2)
        if dist <= 2.5:
            sql = 'insert into comparative_tables_interpro_profiles_euclidian_distance_%s values ("%s", "%s", %s);' % (one_pair[0],
                                                                                                                       one_pair[1],
                                                                                                                       dist)
            server.adaptor.execute(sql,)
            server.adaptor.commit()


def euclidian_dist_interpro(biodb):

    from chlamdb.biosqldb import manipulate_biosqldb
    from multiprocessing import Process
    import math

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_profiles_table = 'CREATE TABLE IF NOT EXISTS comparative_tables_interpro_profiles_euclidian_distance (interpro_1 varchar(100), ' \
                         ' interpro_2 varchar(100), euclidian_dist FLOAT, INDEX interpro_1 (interpro_1), INDEX interpro_2 (interpro_2))'
    try:
        print (sql_profiles_table)
        server.adaptor.execute(sql_profiles_table)
    except:
        print ('problem creating the sql table')

    sql = 'select * from comparative_tables_interpro'
    interpro2profile = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    combinations = []

    # get a list of all combination of groups
    # allows then to split the job into 8 process
    total = len(interpro2profile.keys())
    for i, orthogroup_1 in enumerate(interpro2profile.keys()):
        #if i == 1000:
        #    break
        print ("%s/%s" % (i, total))
        if sum(interpro2profile[orthogroup_1]) == 1:
            print ('skip!')
            continue
        count = 0
        for n in interpro2profile[orthogroup_1]:
            if n>0:
                count+=1
        if count>1:
            print ('range', i, total)
            for y, orthogroup_2 in enumerate(interpro2profile.keys()[i:]):# in range(i, total):
                #orthogroup_2 = orthogroup2profile.keys()[y]
                if sum(interpro2profile[orthogroup_2]) == 1:
                    continue
                count = 0
                for n in interpro2profile[orthogroup_2]:
                    if n>0:
                        count+=1
                if count<=1:
                    continue
                else:
                    combinations.append((orthogroup_1, orthogroup_2))
    n_cpu = 8
    print ('n of pairs:', len(combinations))
    #import sys
    #sys.exit()
    n_poc_per_list = math.ceil(len(combinations)/float(n_cpu))
    #print "n_poc_per_list", n_poc_per_list
    #print range(0, len(combinations))
    # split the jobs into 8 lists
    query_lists = chunks(range(0, len(combinations)), int(n_poc_per_list))
    # start the 8 processes
    procs = []
    for one_list in query_lists:
        comb_list = [combinations[i] for i in one_list]
        #print len(comb_list)
        proc = Process(target=sql_euclidian_dist_interpro, args=(biodb, comb_list, interpro2profile))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-c", '--cpus', type=int, help="Number of cpus", default=8)
    parser.add_argument("-j", '--jaccard', help="Jaccard Index", action="store_true")
    parser.add_argument("-e", '--euclidean', help="Euclidean distance", action="store_true")

    args = parser.parse_args()

    if args.euclidean:
        euclidian_dist_orthogroups(args.db_name, 
                                   args.cpus)
    if args.jaccard:
        jaccard_dist_orthogroups(args.db_name, 
                                args.cpus)
    
    #euclidian_dist_cogs(args.db_name)
    #euclidian_dist_interpro(args.db_name)
