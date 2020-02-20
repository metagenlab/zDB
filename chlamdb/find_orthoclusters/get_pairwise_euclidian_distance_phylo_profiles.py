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
        print ("%s/%s" % (i,n))
        dist = euclidean(orthogroup2profile[one_pair[0]], orthogroup2profile[one_pair[1]])
        if dist <= 2.5:
            sql = 'insert into comparative_tables.phylo_profiles_eucl_dist2_%s values ("%s", "%s", %s);' % (biodb,
                                                                                                        one_pair[0],
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


def euclidian_dist_orthogroups(biodb, merge_taxons=False):

    from chlamdb.biosqldb import manipulate_biosqldb
    from multiprocessing import Process
    import math
    import pandas
    import numpy
    import os
    import rpy2.robjects.numpy2ri
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri as numpy2ri

    rpy2.robjects.numpy2ri.activate()

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_profiles_table = 'CREATE TABLE IF NOT EXISTS comparative_tables.phylo_profiles_eucl_dist2_%s (group_1 varchar(100), ' \
                         ' group_2 varchar(100), euclidian_dist FLOAT, INDEX group_1 (group_1), INDEX group_2 (group_2))' % (biodb)
    try:
        print (sql_profiles_table)
        server.adaptor.execute(sql_profiles_table)
    except:
        print ('problem creating the sql table')

    sql = 'select * from comparative_tables_orthology' % biodb

    sql2 = 'show columns from comparative_tables_orthology' % biodb


    # get matrix as pantas table: orthogroups as rows, genomes as columns
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

    sql = 'select biodatabase_id from biodatabase where name="%s"' % biodb

    biodb_id = server.adaptor.execute_and_fetchall(sql,)[0][0]


    sqlpsw = os.environ['SQLPSW']

    robjects.r("""
    library("RMySQL")
    library(reshape2)

    con <- dbConnect(MySQL(),
             user="root", password="%s",
             dbname="comparative_tables", host="localhost")

    rs1 <- dbSendQuery(con, 'select taxon_1,taxon_2, median_identity from shared_og_av_id_%s union select taxon_2, taxon_1, median_identity from shared_og_av_id_%s;')
    pairwise_identity<- dbFetch(rs1, n=-1)
    rs2 <- dbSendQuery(con, 'select taxon_id,description from biosqldb.bioentry where biodatabase_id=%s and description not like "%%plasmid%%";')
    taxon2description<- dbFetch(rs2, n=-1)

    pairwise_identity_matrix <- dcast(pairwise_identity, taxon_1~taxon_2)
    print(dim(pairwise_identity))
    rownames(pairwise_identity_matrix) <- pairwise_identity_matrix$taxon_1

    pairwise_identity_matrix<-pairwise_identity_matrix[,2:length(pairwise_identity_matrix)]

    pairwise_dist <- as.dist(100-pairwise_identity_matrix)
    print(pairwise_dist)
    hc <- hclust(pairwise_dist)
    clusterCut <- cutree(hc,h=30)
    taxons <- names(clusterCut)

    """ % (sqlpsw, biodb, biodb, biodb_id))

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
    merged_profiles = merge_dataframe_columns(profiles, cluster2taxons)
    # replace multicounts by 1
    merged_profiles[merged_profiles >1] = 1
    #print groups[200]
    combinations = []
    #print merged_profiles.loc[200,:]
    #import sys
    #sys.exit()

    orthogroup2profile = {}
    for i, group in enumerate(groups):
        orthogroup2profile[group] = merged_profiles.loc[i,]

    filtered_groups = []
    total = len(groups)
    for i, orthogroup in enumerate(groups):
        print ("%s/%s" % (i, total))
        #if no homologs in other genomes, skip
        if sum(merged_profiles.loc[i,:]) == 1:
            continue
        else:
            filtered_groups.append(orthogroup)

    # get a list of all combination of groups
    # allows then to split the job into 8 process
    total = len(groups)
    for i, orthogroup_1 in enumerate(filtered_groups):
        print ("%s/%s" % (i, total))
        for y, orthogroup_2 in enumerate(filtered_groups[i:]):# in range(i, total):
            # if no homolog in other genomes, skip
            combinations.append((orthogroup_1, orthogroup_2))
    n_cpu = 8
    # n of pairs: 139259021
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
        proc = Process(target=sql_euclidian_dist_orthogroups, args=(biodb, comb_list, orthogroup2profile))
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
            sql = 'insert into comparative_tables.cogs_profiles_euclidian_distance_%s values ("%s", "%s", %s);' % (biodb,
                                                                                                        one_pair[0],
                                                                                                        one_pair[1],
                                                                                                        dist)
            server.adaptor.execute(sql,)
            server.adaptor.commit()


def euclidian_dist_cogs(biodb):

    from chlamdb.biosqldb import manipulate_biosqldb
    from multiprocessing import Process
    import math

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_profiles_table = 'CREATE TABLE IF NOT EXISTS comparative_tables.cogs_profiles_euclidian_distance_%s (cog_1 varchar(100), ' \
                         ' cog_2 varchar(100), euclidian_dist FLOAT, INDEX cog_1 (cog_1), INDEX cog_2 (cog_2))' % (biodb)
    try:
        print (sql_profiles_table)
        server.adaptor.execute(sql_profiles_table)
    except:
        print ('problem creating the sql table')

    sql = 'select * from comparative_tables_COG' % biodb
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
        print ("%s/%s" % (i,n))
        list1 = np.asarray(interpro2profile[one_pair[0]])
        list2 = np.asarray(interpro2profile[one_pair[1]])
        list1[list1 > 1] = 1
        list2[list2 > 1] = 1

        dist = euclidean(list1, list2)
        if dist <= 2.5:
            sql = 'insert into comparative_tables.interpro_profiles_euclidian_distance_%s values ("%s", "%s", %s);' % (biodb,
                                                                                                        one_pair[0],
                                                                                                        one_pair[1],
                                                                                                        dist)
            server.adaptor.execute(sql,)
            server.adaptor.commit()


def euclidian_dist_interpro(biodb):

    from chlamdb.biosqldb import manipulate_biosqldb
    from multiprocessing import Process
    import math

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_profiles_table = 'CREATE TABLE IF NOT EXISTS comparative_tables.interpro_profiles_euclidian_distance_%s (interpro_1 varchar(100), ' \
                         ' interpro_2 varchar(100), euclidian_dist FLOAT, INDEX interpro_1 (interpro_1), INDEX interpro_2 (interpro_2))' % (biodb)
    try:
        print (sql_profiles_table)
        server.adaptor.execute(sql_profiles_table)
    except:
        print ('problem creating the sql table')

    sql = 'select * from comparative_tables_interpro' % biodb
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

    args = parser.parse_args()

    euclidian_dist_orthogroups(args.db_name)
    #euclidian_dist_cogs(args.db_name)
    #euclidian_dist_interpro(args.db_name)
