#!/usr/bin/python

from BioSQL import BioSeqDatabase
import sys
from Bio import AlignIO
import numpy as np
import numpy
from multiprocessing import Process, Queue
import time
from multiprocessing import cpu_count
import os
import MySQLdb
from chlamdb.biosqldb import manipulate_biosqldb


def check_identity(cursor, 
                   orthogroup, 
                   locus1, 
                   locus2):

    sql = 'select identity from orthology_identity where (locus_a = "%s" and locus_b = "%s") or (locus_a ="%s" and locus_b="%s")'  % (locus1, 
                                                                                                                                      locus2, 
                                                                                                                                      locus1, 
                                                                                                                                      locus2)
    print(sql)
    cursor.execute(sql)
    identity = int(cursor.fetchone()[0])
    return identity


def heatmap_presence_absence(biodb_name, group_name):
    server, db = manipulate_biosqldb.load_db(biodb_name)

    genomes = manipulate_biosqldb.get_genome_taxons_list(server, biodb_name)

    template = ''
    for i in range(0, len(genomes)):
        template += '`%s`, ' % genomes[i]
    template += '`%s`' % genomes[-1]

    #print "template", template

    sql = 'select %s from orthology_%s where orthogroup = "%s"' % (template, biodb_name, group_name)

    result = [int(i) for i in server.adaptor.execute_and_fetchall(sql,)[0]]
    #print result
    taxon2presence_absence = {}
    for x, y in zip(genomes, result):
        taxon2presence_absence[x] = y
    return taxon2presence_absence


def orthogroup2identity_dico(biodb_name, orthogroup):
    server, db = manipulate_biosqldb.load_db(biodb_name)
    sql = 'select locus_a, locus_b, identity from orthology_identity where orthogroup="%s"' % orthogroup
    identity_table = server.adaptor.execute_and_fetchall(sql,)
    locus2locus_id = {}
    for row in identity_table:
        if row[0] not in locus2locus_id:
            locus2locus_id[row[0]] = {}
            locus2locus_id[row[0]][row[1]] = row[2]
        else:
            locus2locus_id[row[0]][row[1]] = row[2]
        if row[1] not in locus2locus_id:
            locus2locus_id[row[1]] = {}
            locus2locus_id[row[1]][row[0]] = row[2]
        else:
            locus2locus_id[row[1]][row[0]] = row[2]
    return locus2locus_id


def locus_list2presence_absence_all_genomes(locus_list, biodb_name):
    server, db = manipulate_biosqldb.load_db(biodb_name)

    locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb_name)

    taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb_name)

    import re
    for i in taxon_id2description.keys():
        taxon_id2description[i] = re.sub(" subsp\. aureus", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(", complete genome\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(", complete sequence\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub("strain ", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub("str\. ", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(" complete genome sequence\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(" complete genome\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(" chromosome", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub("Staphylococcus aureus ", "", taxon_id2description[i])


    header = 'orthogroup\t'
    genomes = manipulate_biosqldb.get_genome_taxons_list(server, biodb_name)
    for i in genomes:
        header += taxon_id2description[i] + '\t'
    final_out = header + '\n'

    for i in locus_list:
        #print "locus", i
        seqfeature_id = locus_tag2seqfeature_id[i]
        orthogroup = manipulate_biosqldb.seqfeature_id2orthogroup(server, seqfeature_id, biodb_name)
        #print "ortho", orthogroup
        dico = heatmap_presence_absence(biodb_name, orthogroup)

        #print "dico done..."
        #print dico
        out = '%s\t' % orthogroup
        for i in genomes:

            out += '%s\t' % dico[i]
        final_out +=out + '\n'

    return final_out



def orthogroup2average_identity(biodatabase_name):
    server, db = manipulate_biosqldb.load_db(biodatabase_name)
    sql = 'select orthogroup,average_identity from orth_%s.average_identity' % biodatabase_name
    result = server.adaptor.execute_and_fetchall(sql)
    return manipulate_biosqldb.to_dict(result)





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-a",'--align_files', type=str, help="aliphment files", nargs='+')
    parser.add_argument("-d",'--db_name', type=str, help="database_name")

    args = parser.parse_args()

    tata = Orthogroup_Identity_DB(args.db_name)
    #print "importing alignments..."
    tata.import_alignments(tata.cursor, args.align_files)
    tata.add_average_orthogroup_identity(args.db_name)
    #check_identity("Chlamydia_12_14", "group_825", "Cav1_00733", "CT565")
