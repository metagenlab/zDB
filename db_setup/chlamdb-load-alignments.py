#!/usr/bin/env python

import sys
from Bio import AlignIO
import numpy as np
import numpy
from multiprocessing import Process, Queue
import time
from multiprocessing import cpu_count
import os

from chlamdb.biosqldb import manipulate_biosqldb
import pickle
import time

class Orthogroup_Identity_DB:
    def __init__(self, database_name, ncpus):

        self.n_cpus = ncpus

        from chlamdb.biosqldb import manipulate_biosqldb

        self.server, self.db = manipulate_biosqldb.load_db(database_name)
        self.conn = self.server.adaptor.conn
        self.cursor = self.server.adaptor.cursor

        self.count = 0


        sql = 'CREATE TABLE if not exists orthology_identity (' \
                   ' id INTEGER PRIMARY KEY,' \
                   ' orthogroup varchar(200),' \
                   ' locus_a VARCHAR(30) NOT NULL,' \
                   ' locus_b VARCHAR(30) NOT NULL,' \
                   ' identity FLOAT(5) )'

        self.cursor.execute(sql)

        sql = f'select locus_tag,orthogroup_name from orthology_seqfeature_id2orthogroup t1 ' \
              f' inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id ' \
              f' inner join orthology_orthogroup t3 on t1.orthogroup_id=t3.orthogroup_id'
              
        self.cursor.execute(sql,)
        self.locus_tag2orthogroup = {}
        for row in self.cursor.fetchall():
            self.locus_tag2orthogroup[row[0]] = row[1]

        sql = 'CREATE TABLE if not exists orthology_average_identity (orthogroup varchar(100), identity float(5))'
        self.server.adaptor.execute(sql)       
        self.server.commit()

    def import_alignments(self, cursor, alignment_files):
        import sys
        sys.stdout.write("getting orthogroup identity matrix from %s fasta alignments...\n" % (len(alignment_files)))
        all_matrix = self._get_group_id2identity_matrix(alignment_files, self.n_cpus)
        sys.stdout.write("Creating mysql table to store identity %s matrix...\n" % len(all_matrix))
        for one_matrix in all_matrix.keys():
            #print one_matrix
            #self._create_identity_table(cursor, all_matrix[one_matrix], one_matrix)
            self._add_identity_data(cursor, all_matrix[one_matrix], one_matrix)
            self.conn.commit()


    def index_identity_table(self):
        
        sql1 = 'create index oidla on orthology_identity(locus_a)'
        sql2 = 'create index oidlb on orthology_identity(locus_b)'
        sql3 = 'create index oidlg on orthology_identity(orthogroup)'
        
        self.cursor.execute(sql1,)
        self.cursor.execute(sql2,)
        self.cursor.execute(sql3,)
        

    def _add_identity_data(self, server, group_matrix, group_name):
        '''
        # old approach creating tables of the size n locus vs n locus
        # not adequate for big matrices: sql tables have a max number of columns
        # ==> use two columns tables: locus 1, locus 1, identity
        '''

        locus_list = [i.decode("utf-8") for i in group_matrix[1:,0]]

        for x in range(1, len(group_matrix[:,0])):
            for y in range(x, len(group_matrix[:,0])):
                if group_matrix[x,y] != 0:
                    sql = 'INSERT INTO orthology_identity (orthogroup, locus_a, locus_b, identity) VALUES ("%s", "%s", "%s", %s)'  % (group_name,
                                                                                                                                      locus_list[x-1],
                                                                                                                                      locus_list[y-1],
                                                                                                                                      float(group_matrix[x,y]))
                else:
                    sql = 'INSERT INTO orthology_identity (orthogroup, locus_a, locus_b, identity) VALUES ("%s", "%s", "%s", NULL)'  % (group_name,
                                                                                                                     locus_list[x-1],
                                                                                                                     locus_list[y-1]
                                                                                                                     )
                self.cursor.execute(sql)

    def _chunks(self, l, n):
        import random
        random.shuffle(l)
        return [l[i:i+n] for i in range(0, len(l), n)]

    def _pairewise_identity(self, seq1, seq2):
        import re
        A=list(seq1)
        B=list(seq2)
        identical_sites = 0
        aligned_sites = 0
        gaps=0
        for n in range(0, len(A)):
            if A[n] != "-" and B[n] != "-":
                aligned_sites+=1
            else:
                continue
            if A[n]==B[n]:
                identical_sites+=1

        # update 08.16
        # calculate coverage
        # chose to cut coverage at 30%
        seq1_no_gap = re.sub('-','', str(seq1.seq))
        seq2_no_gap = re.sub('-','', str(seq2.seq))

        if float(aligned_sites)/len(seq1_no_gap) < 0.3:
            return 0
        elif float(aligned_sites)/len(seq2_no_gap) < 0.3:
            return 0
        else:
            return 100*(identical_sites/float(aligned_sites))


    def _get_identity_matrix_from_multiple_alignment(self, alignment):
        identity_matrix = np.chararray((len(alignment)+1, len(alignment)+1), itemsize=30)
        identity_matrix[0, 0] = "-"
        for x in range(0, len(alignment)):
            identity_matrix[x+1, 0] = alignment[x].name
            identity_matrix[0, x+1] = alignment[x].name

            for y in range(x, len(alignment)):
                identity = self._pairewise_identity(alignment[x], alignment[y])
                identity_matrix[y+1, x+1] = round(identity, 2)
                identity_matrix[x+1, y+1] = round(identity, 2)
        return identity_matrix


    def _group_id(self, group_align_files, out_q, list_id):
        outdict = {}
        for i, align_file in enumerate(group_align_files):
            align = AlignIO.read(align_file, "fasta")
            id_matrix = self._get_identity_matrix_from_multiple_alignment(align)
            group_name = self.locus_tag2orthogroup[align[0].id] #os.path.basename(align_file).split(".")[0]
            outdict[group_name] = id_matrix
            self.count += 1
            if i % 20 == 0:
                print("List %s: %s / %s" % (list_id, i, len(group_align_files)))
        pickle.dump(outdict, open("list_%s.p" % list_id, "wb"))
        time.sleep(2)
        out_q.put("list_%s.p" % list_id)

    def _get_group_id2identity_matrix(self, alignments, n_cpus):

        out_q = Queue(n_cpus)

        n_cpu = n_cpus
        n_poc_per_list = int(numpy.ceil(len(alignments)/float(n_cpu)))
        query_lists = self._chunks(alignments, n_poc_per_list)
        #print query_lists
        procs = []
        print("starting... %s prallel jobs" % n_cpu)
        for n, one_list in enumerate(query_lists):
            print("list", n, len(one_list), "elements")
            proc = Process(target=self._group_id, args=(one_list, out_q, n))
            procs.append(proc)
            proc.start()

        # Collect all results into a single result dict. We know how many dicts
        # with results to expect.
        group_id2identity_matrix = {}
        for i in range(n_cpu):
            group_id2identity_matrix.update(pickle.load(open(out_q.get(),"rb")))
        #print "join proc"
        #print "n groups:", len(group_id2identity_matrix.keys())
        time.sleep(5)


        # Wait for all worker processes to finish
        for proc in procs:
            proc.join()

        return group_id2identity_matrix


    def _get_average_identity_from_identity_matrix(self, id_matrix):
        # get upper part of the matrix indexes
        indexes = np.triu_indices(len(id_matrix))
        # calculate mean
        # remove "False" elements (no alignments between some pairs of the orthogroup
        mask = np.invert(np.in1d(id_matrix[indexes], [False]))
        all_id = id_matrix[indexes][mask]
        return np.mean(all_id)


    def add_average_orthogroup_identity(self):

        sql = f'select orthogroup from comparative_tables_orthology'
                   
        groups = [i[0] for i in self.server.adaptor.execute_and_fetchall(sql, )]
        #print len(groups)

        for group in groups:
            id_values = self.get_orthogroup_identity_table(group)
            if len(id_values) > 0:
                av_id = round(np.mean(id_values), 2)
                
                sql = 'insert into orthology_average_identity values ("%s", %s)' % (group, av_id)

                self.server.adaptor.execute(sql)
        self.server.commit()


    def get_orthogroup_identity_table(self, orthogroup):
        import os

        sql = 'SELECT identity FROM orthology_identity where orthogroup="%s"' % orthogroup

        values = [float(i[0]) for i in self.server.adaptor.execute_and_fetchall(sql,)]
        return values



def locus_tag2identity_best_hit_all_genomes(biodb_name, locus, group_name, locus_tag2taxonomic_id_dict = ""):
    server, db = manipulate_biosqldb.load_db(biodb_name)
    genomes = manipulate_biosqldb.get_genome_taxons_list(server, biodb_name)

    # create empty dictionnary
    genome2genome = {}
    for genome1 in genomes:
        if not genome1 in genome2genome.keys():
            genome2genome[str(genome1)] = {}
            for genome2 in genomes:
                genome2genome[str(genome1)][str(genome2)] = []

    try:
        identity_table = np.array(get_orthogroup_identity_table(biodb_name, group_name))
        all_locus = identity_table[:, 0]
    except:
        genome2best_hit = {}
        for genome1 in genomes:
            genome2best_hit[genome1] = 0
        return genome2best_hit

    locus_tag2taxonomic_id_dict = manipulate_biosqldb.locus_tag2genome_taxon_id(server, biodb_name)
    locus_index_ref = list(all_locus).index(locus)

    '''
    for y in range(1, len(all_locus) + 1):
        genome1 = locus_tag2taxonomic_id_dict[all_locus[y-1]]
        for z in range(0, len(all_locus)):
            genome2= locus_tag2taxonomic_id_dict[all_locus[z]]
            genome2genome[str(genome1)][str(genome2)].append(identity_table[z, y])
    '''

    genome2best_hit = {}
    for genome1 in genomes:
        genome2best_hit[genome1] = 0
    for target_locus in all_locus:
        locus_index_target = list(all_locus).index(target_locus)
        identity = identity_table[locus_index_ref, locus_index_target + 1]
        target_taxon = locus_tag2taxonomic_id_dict[target_locus]
        if float(identity) > genome2best_hit[str(target_taxon)]:
            genome2best_hit[str(target_taxon)] = float(identity)
    return genome2best_hit


def locus_list2identity_in_other_genomes(locus_list, biodb):
    server, db = manipulate_biosqldb.load_db(biodb)
    locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)
    taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

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
    dico = locus_tag2identity_best_hit_all_genomes(biodb, 'wcw_1594', 'group_417')
    for i in dico.keys():

        header += taxon_id2description[i] + '\t'

    final_out = header + '\n'

    for locus in locus_list:
        #print "locus", i
        seqfeature_id = locus_tag2seqfeature_id[locus]
        orthogroup = manipulate_biosqldb.seqfeature_id2orthogroup(server, seqfeature_id, biodb)
        print ("ortho", orthogroup)
        dico = locus_tag2identity_best_hit_all_genomes(biodb, locus, orthogroup)
        print ("dico done...")
        out = '%s\t' % orthogroup
        for i in dico.keys():
            identity = dico[i]
            out += '%s\t' % identity
        final_out += out + '\n'
    return final_out


def heatmap_presence_absence(biodb_name, group_name):
    server, db = manipulate_biosqldb.load_db(biodb_name)
    genomes = manipulate_biosqldb.get_genome_taxons_list(server, biodb_name)
    template = ''
    for i in range(0, len(genomes)):
        template += '`%s`, ' % genomes[i]
    template += '`%s`' % genomes[-1]

    #print "template", template

    sql = f'select %s from orthology where orthogroup = "%s"' % (template, group_name)

    result = [int(i) for i in self.server.adaptor.execute_and_fetchall(sql,)[0]]
    #print result
    taxon2presence_absence = {}
    for x, y in zip(genomes, result):
        taxon2presence_absence[x] = y
    return taxon2presence_absence


def orthogroup2identity_dico(biodb_name, orthogroup):
    identity_table = np.array(get_orthogroup_identity_table(biodb_name, orthogroup))
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
        print ("ortho", orthogroup)
        dico = heatmap_presence_absence(biodb_name, orthogroup)

        out = '%s\t' % orthogroup
        for i in genomes:
            out += '%s\t' % dico[i]
        final_out +=out + '\n'

    return final_out


def orthogroup2average_identity(biodatabase_name):
    server, db = manipulate_biosqldb.load_db(biodatabase_name)
    sql = 'select orthogroup, average_identity from orthology_average_identity'
    result = self.server.adaptor.execute_and_fetchall(sql)
    return manipulate_biosqldb.to_dict(result)


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-a",'--align_files', type=str, help="aliphment files", nargs='+')
    parser.add_argument("-d",'--db_name', type=str, help="database_name")
    parser.add_argument("-c",'--cpus', type=int, help="n cpus")

    args = parser.parse_args()

    tata = Orthogroup_Identity_DB(args.db_name, args.cpus)
    print("importing alignments...")
    
    tata.import_alignments(tata.cursor, args.align_files)  
    print("Index table")
    tata.index_identity_table()
    print("add_average_orthogroup_identity")
    tata.add_average_orthogroup_identity()
    
    manipulate_biosqldb.update_config_table(args.db_name, "orthogroup_alignments")
    
