#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def create_DOOR_operon_table(biodb):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    # OperonID	GI	Synonym	Start	End	Strand	Length	COG_number	Product
    sql = 'CREATE TABLE IF NOT EXISTS custom_tables_DOOR2_operons (operon_id INT, gi INT, seqfeature_id INT, old_locus_tag varchar(400), COG_number ' \
          ' varchar(400), product TEXT, index seqfeature_id (seqfeature_id), index old_locus_tag(old_locus_tag), index operon_id(operon_id))'
    server.adaptor.execute(sql)
    server.commit()

def door_accession2door_operon_table(accession):

    import urllib2

    door_link = 'http://csbl.bmb.uga.edu/DOOR/downloadNCoperon.php?NC_id=%s' % accession

    req = urllib2.Request(door_link)
    data = urllib2.urlopen(req)

    operon_table = []
    try:
        for i, row in enumerate(data):
            if i != 0 :
                data = row.rstrip().split('\t')
                if len(data) == 9:
                    operon_table.append(data)
        return operon_table
    except:
        print ('echec')
        import time
        time.sleep(10)
        return door_accession2door_operon_table(accession)


def accession2operon_table(biodb):

    from chlamdb.biosqldb import manipulate_biosqldb
    import re

    server, db = manipulate_biosqldb.load_db(biodb)

    create_DOOR_operon_table(biodb)

    accession_list_sql = 'select t1.accession from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                         ' where t2.name="%s"' % biodb
    accession_list = [i[0] for i in server.adaptor.execute_and_fetchall(accession_list_sql,)]

    # get full join
    sql = 'select old_locus_tag,t1.seqfeature_id from custom_tables_locus2seqfeature_id t1' \
          ' inner join custom_tables_seqfeature_id2old_locus_tag t2 ' \
          ' on t1.seqfeature_id=t2.seqfeature_id;'
    print (sql)

    try:
        old_locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    except:
        old_locus2seqfeature_id = {}
    print ("old_locus2seqfeature_id", old_locus2seqfeature_id)
    sql = 'select locus_tag, seqfeature_id from custom_tables.locus2seqfeature_id'

    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    #print locus_tag2seqfeature_id

    pattern_ctracho = 'CT([0-9]+)'

    for accession in accession_list:
        print("accession", accession)
        door_id_sql = 'select door_id from custom_tables_DOOR2_operon_accessions where accession="%s"' % accession

        try:
            door_id = server.adaptor.execute_and_fetchall(door_id_sql,)[0][0]
        except IndexError:
            print('%s Not in doord database, skipping...' % accession)
            continue
        print (door_id)
        door_data = door_accession2door_operon_table(door_id)
        for door_entry in door_data:

            m = re.match(pattern_ctracho, door_entry[2])

            # C trachomatis special case: add underscore

            if m:
                door_entry[2] = 'CT_%s' % m.group(1)
            try:
                new_locus_seqfeature_id = old_locus2seqfeature_id[door_entry[2]]
            except KeyError:
                try:
                    new_locus_seqfeature_id = locus_tag2seqfeature_id[door_entry[2]]
                except KeyError:
                    print (door_entry[2])
                    new_locus_seqfeature_id = 0
            sql = 'INSERT INTO custom_tables_DOOR2_operons values (%s, %s, %s, "%s", "%s", "%s")' % (door_entry[0],
                                                                                                     door_entry[1],
                                                                                                     new_locus_seqfeature_id,
                                                                                                     door_entry[2],
                                                                                                     door_entry[7],
                                                                                                     door_entry[8])

            server.adaptor.execute(sql,)
            server.commit()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--biodb', type=str, help="biodatabase name")
    parser.add_argument("-c", '--create_accession_table', action='store_true', help="biodatabase name")

    args = parser.parse_args()
    if args.create_accession_table:
        door_accession2door_operon_table(accession=1101)
    accession2operon_table(args.biodb)
