#! /usr/bin/env python


import sys
from Bio import GenBank
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from manipulate_biosqldb import load_db
from manipulate_biosqldb import query_yes_no


def chunks(l, n):
    "return sublists of l of minimum length n (work subdivision for the subprocesing module"
    return [l[i:i+n] for i in range(0, len(l), n)]



def sql_euclidian_dist_orthogroups(biodb, one_list, orthogroup2profile):
    from scipy.spatial.distance import euclidean
    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    n = len(one_list)
    for i, one_pair in enumerate(one_list):
        print "%s/%s" % (i,n)
        dist = euclidean(orthogroup2profile[one_pair[0]], orthogroup2profile[one_pair[1]])
        if dist <= 2.5:
            sql = 'insert into comparative_tables.phylogenetic_profiles_euclidian_distance_%s values ("%s", "%s", %s);' % (biodb,
                                                                                                        one_pair[0],
                                                                                                        one_pair[1],
                                                                                                        dist)
            server.adaptor.execute(sql,)
            server.adaptor.commit()


def euclidian_dist_orthogroups(biodb):

    import manipulate_biosqldb
    from multiprocessing import Process
    import math

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_profiles_table = 'CREATE TABLE IF NOT EXISTS comparative_tables.phylogenetic_profiles_euclidian_distance_%s (group_1 varchar(100), ' \
                         ' group_2 varchar(100), euclidian_dist FLOAT, INDEX group_1 (group_1), INDEX group_2 (group_2))' % (biodb)
    try:
        print sql_profiles_table
        server.adaptor.execute(sql_profiles_table)
    except:
        print 'problem creating the sql table'

    sql = 'select * from orthology_%s' % biodb
    orthogroup2profile = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    combinations = []

    # get a list of all combination of groups
    # allows then to split the job into 8 process
    total = len(orthogroup2profile.keys())
    for i, orthogroup_1 in enumerate(orthogroup2profile.keys()):
        #if i == 1000:
        #    break
        print "%s/%s" % (i, total)
        if sum(orthogroup2profile[orthogroup_1]) == 1:
            print 'skip!'
            continue
        count = 0
        for n in orthogroup2profile[orthogroup_1]:
            if n>0:
                count+=1
        if count>1:
            print 'range', i, total
            for y, orthogroup_2 in enumerate(orthogroup2profile.keys()[i:]):# in range(i, total):
                #orthogroup_2 = orthogroup2profile.keys()[y]
                if sum(orthogroup2profile[orthogroup_2]) == 1:
                    continue
                count = 0
                for n in orthogroup2profile[orthogroup_2]:
                    if n>0:
                        count+=1
                if count<=1:
                    continue
                else:
                    combinations.append((orthogroup_1, orthogroup_2))
    n_cpu = 8
    print 'n of pairs:', len(combinations)
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
    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    n = len(one_list)
    for i, one_pair in enumerate(one_list):
        print "%s/%s" % (i,n)
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

    import manipulate_biosqldb
    from multiprocessing import Process
    import math

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_profiles_table = 'CREATE TABLE IF NOT EXISTS comparative_tables.cogs_profiles_euclidian_distance_%s (cog_1 varchar(100), ' \
                         ' cog_2 varchar(100), euclidian_dist FLOAT, INDEX cog_1 (cog_1), INDEX cog_2 (cog_2))' % (biodb)
    try:
        print sql_profiles_table
        server.adaptor.execute(sql_profiles_table)
    except:
        print 'problem creating the sql table'

    sql = 'select * from comparative_tables.COG_%s' % biodb
    cog2profile = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    combinations = []

    # get a list of all combination of groups
    # allows then to split the job into 8 process
    total = len(cog2profile.keys())
    for i, orthogroup_1 in enumerate(cog2profile.keys()):
        #if i == 1000:
        #    break
        print "%s/%s" % (i, total)
        if sum(cog2profile[orthogroup_1]) == 1:
            print 'skip!'
            continue
        count = 0
        for n in cog2profile[orthogroup_1]:
            if n>0:
                count+=1
        if count>1:
            print 'range', i, total
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
    print 'n of pairs:', len(combinations)
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
    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    n = len(one_list)
    for i, one_pair in enumerate(one_list):
        print "%s/%s" % (i,n)
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

    import manipulate_biosqldb
    from multiprocessing import Process
    import math

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_profiles_table = 'CREATE TABLE IF NOT EXISTS comparative_tables.interpro_profiles_euclidian_distance_%s (interpro_1 varchar(100), ' \
                         ' interpro_2 varchar(100), euclidian_dist FLOAT, INDEX interpro_1 (interpro_1), INDEX interpro_2 (interpro_2))' % (biodb)
    try:
        print sql_profiles_table
        server.adaptor.execute(sql_profiles_table)
    except:
        print 'problem creating the sql table'

    sql = 'select * from comparative_tables.interpro_%s' % biodb
    interpro2profile = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    combinations = []

    # get a list of all combination of groups
    # allows then to split the job into 8 process
    total = len(interpro2profile.keys())
    for i, orthogroup_1 in enumerate(interpro2profile.keys()):
        #if i == 1000:
        #    break
        print "%s/%s" % (i, total)
        if sum(interpro2profile[orthogroup_1]) == 1:
            print 'skip!'
            continue
        count = 0
        for n in interpro2profile[orthogroup_1]:
            if n>0:
                count+=1
        if count>1:
            print 'range', i, total
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
    print 'n of pairs:', len(combinations)
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
    euclidian_dist_cogs(args.db_name)
    euclidian_dist_interpro(args.db_name)
