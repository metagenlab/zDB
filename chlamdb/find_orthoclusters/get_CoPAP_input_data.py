#!/usr/bin/python

def get_profile_fasta(biodb, taxon_id):
    '''

    - ordered taxon
    - transposed orthology table => each row is a different taxon

    :return:
    '''

    import manipulate_biosqldb
    import pandas
    import numpy
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select taxon_id, accession from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where (t1.name="%s" and t2.description not like "%%%%plasmid%%%%")' % biodb


    taxon2accession = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    taxon_id_filter = '`'+'`,`'.join(taxon2accession.keys())+'`'

    sql = 'select t2.locus_tag,%s from comparative_tables_orthology t1 inner join orthology_detail t2 on t1.orthogroup=t2.orthogroup' \
          ' where t2.taxon_id=%s' % (taxon_id_filter, biodb, biodb, taxon_id)
    sql3 = 'show columns from comparative_tables_orthology' % (biodb)

    data = numpy.array([list(i) for i in server.adaptor.execute_and_fetchall(sql,)])

    all_cols = [i[0] for i in server.adaptor.execute_and_fetchall(sql3,)]

    count_df = pandas.DataFrame(data, columns=all_cols)

    count_df = count_df.set_index(['orthogroup'])
    count_df = count_df.apply(pandas.to_numeric, args=('coerce',))
    count_df[(count_df > 1)] = 1
    #print count_df
    transposed_table = count_df.transpose()
    #print transposed_table
    #transposed_table.columns = []
    all_records = []
    for taxon, row in transposed_table.iterrows():
        #print taxon, row
        profile_dat = [str(i) for i in row]
        profile = ''.join(profile_dat)
        simple_seq = Seq(profile)
        simple_seq_r = SeqRecord(simple_seq)
        simple_seq_r.id = taxon2accession[taxon]
        simple_seq_r.description = ""
        all_records.append(simple_seq_r)
    with open("profiles_all.fasta", 'w') as tt:
        SeqIO.write(all_records, tt, 'fasta')
    #print taxon, len(profile)

        #print ">%s\n%s" % (taxon2accession[taxon], profile)

def get_annotation_file(biodb, taxon_id):
    '''
    Pos	ID	Description
    1	G0001	Module A subunit 1
    2	G0002	Module A subunit 2

    :return:
    '''
    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select t2.locus_tag,t2.product from comparative_tables_orthology t1 inner join orthology_detail t2 on t1.orthogroup=t2.orthogroup' \
          ' where t2.taxon_id=%s' % (biodb, biodb, taxon_id)
    data = server.adaptor.execute_and_fetchall(sql,)
    print 'Pos\tID\tDescription'
    for i, row in enumerate(data):
        print "%s\t%s\t%s" % (i+1, row[0], row[1])

def get_phylogeny(biodb):
    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]
    return tree

#get_profile_fasta('chlamydia_04_16', 66)
#print get_phylogeny('chlamydia_04_16')
get_annotation_file('chlamydia_04_16', 66)