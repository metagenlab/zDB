#!/usr/bin/python

from chlamdb.biosqldb import manipulate_biosqldb
from chlamdb.plots import gbk2circos

def locus_tag2orthogroup_size(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)
    sql='select locus_tag, count(*) from orthology_detail group by orthogroup'

    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

def superkingdom_table():
    pass

def get_orthology_matrix_merging_plasmids(server, biodatabase_name, taxon_list=False):


    server, db = manipulate_biosqldb.load_db(biodatabase_name)
    sql='select orthogroup, count(*) from orthology_detail group by orthogroup' % biodatabase_name

    all_orthogroups = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    if not taxon_list:
        all_taxons = manipulate_biosqldb.get_taxon_id_list(server, biodatabase_name)
    else:
        all_taxons = taxon_list

    #print "taxons"
    #print all_taxons
    #print 'get dico ortho count'
    detailed_orthology_count = {}
    for group in all_orthogroups.keys():
        detailed_orthology_count[group] = {}
        for taxon_id in all_taxons:
            detailed_orthology_count[group][int(taxon_id)]= 0


    for taxon_id in all_taxons:
        #print taxon_id
        sql = "select orthogroup, `%s` from comparative_tables_orthology;" % (taxon_id, biodatabase_name)
        #print sql
        dico = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        #print dico
        #import time
        #time.sleep(3)
        #print "taxon", taxon_id
        for group in dico.keys():
            detailed_orthology_count[group][int(taxon_id)] += dico[group]
    #print 'count ok'
    ##print detailed_orthology_count
    return detailed_orthology_count


def COG_tables():
    pass

'''


<domain-id>, <genome-name>, <protein-id>,<protein-length>,
<domain-start>, <domain-end>, <COG-id>, <membership-class>,

* Example:

333894695,Alteromonas_SN2_uid67349,333894695,427,1,427,COG0001,0,

CREATE TABLE cog_2014 (domain_id varchar(100),
                        genome_name varchar(200),
                        protein_id INT,
                        protein_length INT,
                        domain_start INT,
                        domain_end INT,
                        COG_id INT,
                        membership_class INT,
                        index COG_id(COG_id))

# COG	func	name

CREATE table cog_names_2014 (COG_id INT,
                             COG_name varchar(100),
                             function varchar(10),
                             name varchar(200),
                             index COG_id(COG_id))

create table locus_tag2COG_chlamydia_03_15 (locus_tag varchar(100), COG_id varchar (100))


'''


def create_contig_table(db_name):
    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'CREATE TABLE contigs_%s(accession VARCHAR(100) NOT NULL, ' \
          ' contig_name VARCHAR(100) UNIQUE, ' \
          ' start INT, ' \
          ' end INT)' % db_name


    #server.adaptor.execute(sql,)

    sql2 = 'select accession from bioentry inner join biodatabase on biodatabase.biodatabase_id = bioentry.biodatabase_id' \
           ' where biodatabase.name = "%s"' % db_name
    #print sql2
    accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql2,)]



    for accession in accessions:
        record = db.lookup(accession=accession)
        #print record.id
        draft_data = gbk2circos.circos_fasta_draft_misc_features(record)
        if len(draft_data) == 0:
            sql = 'INSERT into contigs_%s (accession, contig_name, start, end) VALUES ("%s", "%s", %s, %s)' % (db_name,
                                                                                                               accession,
                                                                                                                accession,
                                                                                                                0,
                                                                                                                len(record.seq))
            server.adaptor.commit()

            server.adaptor.execute(sql,)
        else:
            for contig in draft_data:
                sql = 'INSERT into contigs_%s (accession, contig_name, start, end) VALUES ("%s", "%s", %s, %s)' % (db_name,
                                                                                                               accession,
                                                                                                               contig[0],
                                                                                                               contig[1],
                                                                                                               contig[2])
                server.adaptor.execute(sql,)
                server.adaptor.commit()


def accession2n_contigs(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select t1.accession, count(*) as n_contigs from bioentry as t1' \
    ' inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id ' \
    ' inner join contigs_chlamydia_03_15 as t3 on t1.accession=t3.accession ' \
    ' where t2.name="%s" group by t3.accession;' % (db_name)

    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


def locus_tag2best_hit(db_name, accession, hit_number=1, rank=False, taxon_name=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    if rank and taxon_name:
        sql_rank_filter = 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.query_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end ' \
                          ' from blastnr_blastnr_hits_%s as t1 ' \
                          ' inner join blastnr_blastnr_hits_taxonomy_filtered_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr_blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr_blastnr_hsps_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.%s = "%s";' % (db_name,
                                                                         accession,
                                                                         db_name,
                                                                         accession,
                                                                         db_name,
                                                                         accession,
                                                                         hit_number,
                                                                         rank,
                                                                         taxon_name)
        #print sql_rank_filter
        data = server.adaptor.execute_and_fetchall(sql_rank_filter,)

        locus_tag2hit = {}

        for one_hsp in data:
            if not one_hsp[0] in locus_tag2hit:
                locus_tag2hit[one_hsp[0]] = [[ i for i in one_hsp[1:]]]
            else:
                locus_tag2hit[one_hsp[0]].append([ i for i in one_hsp[1:]])

        return locus_tag2hit


    else:
        sql_no_rank_filter = 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.query_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end ' \
                          ' from blastnr_blastnr_hits_%s as t1 ' \
                          ' inner join blastnr_blastnr_hits_taxonomy_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr_blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr_blastnr_hsps_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s' % (db_name,
                                                       accession,
                                                       db_name,
                                                       accession,
                                                       db_name,
                                                       accession,
                                                       hit_number)

        #print sql_no_rank_filter
        data = server.adaptor.execute_and_fetchall(sql_no_rank_filter,)

        locus_tag2hit = {}

        for one_hsp in data:
            if not one_hsp[0] in locus_tag2hit:
                locus_tag2hit[one_hsp[0]] = [list(one_hsp[1:])]
            else:
                locus_tag2hit[one_hsp[0]].append(list(one_hsp[1:]))

        return locus_tag2hit

def locus_tag2best_hit_n_taxon_ids(db_name, accession):
    server, db = manipulate_biosqldb.load_db(db_name)
    sql = 'select t3.locus_tag, count(*) as n_taxons from blastnr_blastnr_hits_taxonomy_%s as t1 ' \
          ' inner join blastnr_blastnr_hits_%s as t3 on t1.nr_hit_id=t3.nr_hit_id ' \
          ' where t3.hit_number=1 group by t3.nr_hit_id' % (db_name, accession, db_name, accession)
    #print sql
    all_blastnr_taxons = server.adaptor.execute_and_fetchall(sql,)


def get_locus2plasmid_or_not(biodb):

    sql = 'SELECT locus_tag, IF(organism like "%%%%plasmid%%%%", "True", "False") AS NewResult from orthology_detail;' % biodb
    #print sql
    server, db = manipulate_biosqldb.load_db(biodb)

    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))



def circos_locus2taxon_highest_identity(biodb,
                                        reference_taxon_id,
                                        use_identity_closest_homolog2_table=False):
    '''
    Given one reference taxon, get a dictionnary of the highest identity of homolog(s) in other taxons

    :param reference_accession: reference taxon_id
    :return:
    '''
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    if not use_identity_closest_homolog2_table:
        sql = 'select locus_tag, orthogroup from orthology_detail where taxon_id = %s' % (biodb,
                                                                                             reference_taxon_id)
        sql2 = 'select locus_tag, taxon_id from orthology_detail' % (biodb)

        # get all locus tags and orthogroups from reference genome
        reference_orthogroup2locus_tag = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
        locus2locus_identity = {}
        i = 0
        # for each locus tag, get highest identity in each other taxons
        for locus_A in reference_orthogroup2locus_tag:
            i+=1
            sql = 'select locus_b,identity from orth_%s.%s where locus_a ="%s" ' \
                  ' UNION select locus_a,identity from orth_%s.%s where locus_b ="%s";' % (biodb,
                                                                                          reference_orthogroup2locus_tag[locus_A],
                                                                                          locus_A,
                                                                                          biodb,
                                                                                          reference_orthogroup2locus_tag[locus_A],
                                                                                          locus_A)
            try:
                locus2identity = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
                locus2locus_identity[locus_A] = {}
                for locus_B in locus2identity:
                    if locus_tag2taxon_id[locus_B] not in locus2locus_identity[locus_A]:
                        locus2locus_identity[locus_A][locus_tag2taxon_id[locus_B]] = [locus2identity[locus_B], locus_B]
                    else:
                        # if identity of a second homolog is higher, use this value
                        if locus2identity[locus_B] > locus2locus_identity[locus_A][locus_tag2taxon_id[locus_B]][0]:
                            locus2locus_identity[locus_A][locus_tag2taxon_id[locus_B]] = [locus2identity[locus_B], locus_B]
            except:
                #print 'no homologs for %s' % locus_A
                pass
        #print '------------------', locus2locus_identity.keys()[0], locus2locus_identity[locus2locus_identity.keys()[0]]

    else:
        sql = 'select t2.locus_tag, t3.locus_tag, identity, taxon_2 from comparative_tables_identity_closest_homolog2 t1 ' \
              ' inner join custom_tables_locus2seqfeature_id t2 on t1.locus_1=t2.seqfeature_id ' \
              ' inner join custom_tables_locus2seqfeature_id t3 on t1.locus_2=t3.seqfeature_id ' \
              ' where taxon_1=%s' % (biodb, biodb, biodb, reference_taxon_id)
        data = server.adaptor.execute_and_fetchall(sql,)
        locus2locus_identity = {}
        for row in data:
            if row[0] not in locus2locus_identity:
                locus2locus_identity[row[0]] = {}
                locus2locus_identity[row[0]][row[3]] = [row[2],row[1]]
            else:
                locus2locus_identity[row[0]][row[3]] = [row[2],row[1]]
    #print '------------------', locus2locus_identity.keys()[0], locus2locus_identity[locus2locus_identity.keys()[0]]
    return locus2locus_identity


def taxon_subset2core_orthogroups(biodb, taxon_list, type="nucleotide", mypath="./"):

    from chlamdb.biosqldb import mysqldb_load_mcl_output
    import sys
    from Bio import SeqIO
    import os
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    server, db = manipulate_biosqldb.load_db(biodb)

    sql_include = ''
    if len(taxon_list) > 0:
        for i in range(0, len(taxon_list)-1):
            sql_include += ' `%s` = 1 and ' % taxon_list[i]
        sql_include+='`%s` = 1' % taxon_list[-1]

    sql ='select orthogroup from comparative_tables_orthology where %s' % (biodb, sql_include)
    #print sql
    sys.stdout.write("getting core orthogroup list\n")
    match_groups = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    sys.stdout.write("N single copy core orthogroups: %s\n" % len(match_groups))
    sys.stdout.write("getting locus tag 2 taxon_id\n")
    sql2 = 'select locus_tag, taxon_id from orthology_detail' % (biodb)
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    if type == "nucleotide":

        sys.stdout.write("getting locus tag 2 sequence\n")
        locus_tag2nucl_sequence = mysqldb_load_mcl_output.locus_tag2nucl_sequence_dict(server, db, biodb)


        sys.stdout.write("writing one fasta/group\n")
        for group in match_groups:
            sql = 'select locus_tag from orthology_detail where orthogroup="%s" and taxon_id in (%s)' % (biodb,
                                                                                                                        group,
                                                                                                                        ','.join([str(i) for i in taxon_list]))
            locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            #print group, locus_list
            seqs = []
            for locus in locus_list:
                if locus_tag2taxon_id[locus] in taxon_list:
                    seq = locus_tag2nucl_sequence[locus]
                    seqs.append(seq)
            path = os.path.join(mypath, group + "_nucl.txt")
            #print path
            with open(path, "w") as f:
                SeqIO.write(seqs, f, "fasta")
    else:

        locus_tag2aa_sequence = mysqldb_load_mcl_output.locus_tag2aa_sequence_dict(server, db, biodb)

        #print 'Number of core groups: %s' % len(match_groups)
        for group in match_groups:
            sql = 'select locus_tag from orthology_detail where orthogroup="%s" and taxon_id in (%s)' % (biodb,
                                                                                                           group,
                                                                                                           ','.join([str(i) for i in taxon_list]))
            locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            #print group, locus_list
            seqs = []
            for locus in locus_list:
                if str(locus_tag2taxon_id[locus]) in taxon_list:
                    #print 'OKKKKKKKKKKKKKKKKKKKKKKK'
                    seq = SeqRecord(Seq(locus_tag2aa_sequence[locus],
                                       IUPAC.protein),
                                       id=str(locus_tag2taxon_id[locus]), name="",
                                       description="")
                    seqs.append(seq)
            path = os.path.join(mypath, group + "_aa.txt")
            #print path
            with open(path, "w") as f:
                SeqIO.write(seqs, f, "fasta")


def group2gene(biodb, group_list, rank_limit=2):

    server, db = manipulate_biosqldb.load_db(biodb)
    group_filter = '","'.join(group_list)
    sql = 'select orthogroup_name,description,count from orthology_orthogroup2gene t1 inner join orthology.orthogroup_%s t2 on t1.group_id=t2.orthogroup_id where rank<%s and orthogroup_name in ("%s");' % (biodb, biodb, rank_limit, group_filter)
    data = server.adaptor.execute_and_fetchall(sql,)
    ortho2gene = {}

    for row in data:
        group = row[0]
        gene = row[1]
        count = row[2]
        if group not in ortho2gene:
            ortho2gene[group] = {}
            ortho2gene[group][gene] = {"count": count}
        else:
            ortho2gene[group][gene] = {"count": count}
    return ortho2gene


def group2product(biodb, group_list, rank_limit=2):

    server, db = manipulate_biosqldb.load_db(biodb)
    group_filter = '","'.join(group_list)
    sql = 'select orthogroup_name,description,count from orthology.orthogroup2product_%s t1 inner join orthology.orthogroup_%s t2 on t1.group_id=t2.orthogroup_id where rank<%s and orthogroup_name in ("%s");' % (biodb, biodb, rank_limit, group_filter)
    data = server.adaptor.execute_and_fetchall(sql,)
    ortho2product = {}

    for row in data:
        group = row[0]
        product = row[1]
        count = row[2]
        if group not in ortho2product:
            ortho2product[group] = {}
            ortho2product[group][product] = {"count": count}
        else:
            ortho2product[group][product] = {"count": count}
    return ortho2product


def group2cog(biodb, group_list, rank_limit=2):

    server, db = manipulate_biosqldb.load_db(biodb)
    group_filter = '","'.join(group_list)
    sql = 'select orthogroup_name,COG_name,t3.description,count,code,t5.description from orthology.orthogroup2cog_%s t1 inner join orthology.orthogroup_%s t2 on t1.group_id=t2.orthogroup_id inner join COG_cog_names_2014 t3 on t1.COG_id=t3.COG_id inner join COG_cog_id2cog_category t4 on t1.COG_id=t4.COG_id inner join COG_code2category t5 on t4.category_id=t5.category_id where rank<%s and orthogroup_name in ("%s");' % (biodb, biodb, rank_limit, group_filter)
    data = server.adaptor.execute_and_fetchall(sql,)
    ortho2cog = {}

    for row in data:
        group = row[0]
        cog_name = row[1]
        cog_description = row[2]
        count = row[3]
        category_id = row[4]
        category_description = row[5]
        if group not in ortho2cog:
            ortho2cog[group] = {}
            ortho2cog[group][cog_name] = {"count": count, "cog_description": cog_description, "category_id": category_id, "category_description": category_description}
        else:
            ortho2cog[group][cog_name] = {"count": count, "cog_description": cog_description, "category_id": category_id, "category_description": category_description}
    return ortho2cog


def group2pfam(biodb, group_list, rank_limit=2):

    server, db = manipulate_biosqldb.load_db(biodb)
    group_filter = '","'.join(group_list)
    sql = 'select orthogroup_name,signature_accession,signature_description,count from orthology.orthogroup2pfam_%s t1 inner join orthology.orthogroup_%s t2 on t1.group_id=t2.orthogroup_id inner join interpro_signature t3 on t1.signature_id=t3.signature_id where rank<%s and orthogroup_name in ("%s");' % (biodb, biodb, rank_limit, group_filter)
    data = server.adaptor.execute_and_fetchall(sql,)
    ortho2pfam = {}

    for row in data:
        group = row[0]
        pfam_name = row[1]
        pfam_description = row[2]
        count = row[3]
        if group not in ortho2pfam:
            ortho2pfam[group] = {}
            ortho2pfam[group][pfam_name] = {"count": count, "pfam_description": pfam_description}
        else:
            ortho2pfam[group][pfam_name] = {"count": count, "pfam_description": pfam_description}
    return ortho2pfam


def orthogroup_list2detailed_annotation(ordered_orthogroups, biodb, taxon_filter=False, accessions=False):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    group_filter = '"' + '","'.join(ordered_orthogroups) + '"'


    columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
              'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
    if not taxon_filter:
        sql_2 = 'select %s from orthology_detail where orthogroup in (%s)' % (columns, biodb, group_filter)
    else:
        if not accessions:
            taxon_filter = [str(i) for i in taxon_filter]
            taxon_f = ','.join(taxon_filter)
            sql_2 = 'select %s from orthology_detail where orthogroup in (%s) and taxon_id in (%s)' % (columns,
                                                                                                          biodb,
                                                                                                          group_filter,
                                                                                                          taxon_f)
        else:
            taxon_filter = [str(i) for i in taxon_filter]
            taxon_f = '"'+'","'.join(taxon_filter)+'"'
            sql_2 = 'select %s from orthology_detail where orthogroup in (%s) and accession in (%s)' % (columns,
                                                                                                           biodb,
                                                                                                           group_filter,
                                                                                                           taxon_f)
    raw_data = server.adaptor.execute_and_fetchall(sql_2,)

    orthogroup2genes = group2gene(biodb, ordered_orthogroups)
    orthogroup2products = group2product(biodb, ordered_orthogroups)
    orthogroup2cogs = group2cog(biodb, ordered_orthogroups)
    orthogroup2pfam = group2pfam(biodb, ordered_orthogroups)

    match_groups_data = []
    for i, group in enumerate(ordered_orthogroups):
        genes_data = ''
        for gene in orthogroup2genes[group]:
            genes_data += '%s (%s)<br/>' % (gene, orthogroup2genes[group][gene]["count"])
        genes_data = genes_data[0:-5]
        product_data = ''
        for product in orthogroup2products[group]:
            product_data += '%s (%s)<br/>' % (product, orthogroup2products[group][product]["count"])
        cog_data = ''
        try:
            for cog in orthogroup2cogs[group]:
                cog_category_id = orthogroup2cogs[group][cog]["category_id"]
                cog_category_description = orthogroup2cogs[group][cog]["category_description"]
                cog_data += '<a href=/fam/%s/cog>' \
                            '%s: %s: %s (%s)</a><br/>' % (cog, cog, cog_category_id, cog_category_description, orthogroup2cogs[group][cog]["count"])
            cog_data = cog_data[0:-5]
        except KeyError:
            cog_data += ' <p>-</p> '

        pfam_data = ''
        try:
            for pfam in orthogroup2pfam[group]:
                pfam_data += '<a href=/fam/%s/pfam>%s: %s (%s)</a><br/>' % (pfam,
                                                                                                pfam,
                                                                                                orthogroup2pfam[group][pfam]["pfam_description"],
                                                                                                orthogroup2pfam[group][pfam]["count"])
            pfam_data = pfam_data[0:-5]
        except KeyError:
            pfam_data += ' <p>-</p> '

        match_groups_data.append([i, group, genes_data, product_data, cog_data, pfam_data])
        n = 1
        extract_result = []
        for one_hit in raw_data:
            extract_result.append((n,) + one_hit)
            n += 1
    return match_groups_data, extract_result


def accession2coding_density(biodb, sqlite=False):
    from chlamdb.biosqldb import manipulate_biosqldb
    from Bio.SeqUtils import GC123

    server, db = manipulate_biosqldb.load_db(biodb, sqlite=sqlite)

    # gc content genes
    sql1 = 'select distinct accession from bioentry t1 ' \
           ' inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id' \
           ' where t2.name="%s"' % biodb

    accession_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql1,)]

    '''
    sql_head = 'create table IF NOT EXISTS annotation.gc_content_%s (taxon_id INT, ' \
          ' seqfeature_id INT,' \
          ' seq_length INT, gc_percent FLOAT, gc_1 FLOAT, gc_2 FLOAT, gc_3 FLOAT)' % biodb

    server.adaptor.execute(sql_head, )

    sql1 = 'create index annotation.seqfeature_gcidx_%s ON gc_content_%s (seqfeature_id);' % (biodb, biodb)
    sql2 = 'create index annotation.taxon_id_gcidx%s ON gc_content_%s (taxon_id);' % (biodb, biodb)
    server.adaptor.execute(sql1, )
    server.adaptor.execute(sql2, )
    '''

    accession2density = {}
    accession2n_trna = {}
    accession2n_rrna = {}
    accession2big_contig_length = {}
    accession2n_contigs_without_cds = {}
    accession2n_countigs_without_BBH_chlamydiae = {}
    for n, accession in enumerate(accession_list):
        print(n, accession)

        # get list of BBH chlamydiae
        sql = 'select locus_tag from blastnr_blastnr t1 ' \
              ' inner join biosqldb.bioentry t2 on t1.query_bioentry_id=t2.bioentry_id ' \
              ' inner join biosqldb.biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
              ' inner join blastnr_blastnr_taxonomy t4 on t1.subject_taxid=t4.taxon_id ' \
              ' inner join custom_tables_locus2seqfeature_id t5 ' \
              ' on t1.seqfeature_id=t5.seqfeature_id ' \
              ' where t1.hit_number=1 and t3.name="%s" and t4.phylum="Chlamydiae" and t2.accession="%s";' % (biodb,
                                                                                                             biodb,
                                                                                                             biodb,
                                                                                                             accession)
        try:
            locus_BBH_chlamydiae = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
        except:
            #print 'echec'
            locus_BBH_chlamydiae = []
        #print sql
        #print 'locus_BBH_chlamydiae', len(locus_BBH_chlamydiae)
        accession2n_contigs_without_cds[accession] = 0
        accession2n_countigs_without_BBH_chlamydiae[accession] = 0
        accession = accession.split('.')[0]
        record = db.lookup(accession=accession)
        long_record_list = []
        start = 0
        stop = 0
        big_contig_length = 0
        for feature in record.features:
            if feature.type == 'assembly_gap':
                # assembly gap within scaffolds
                #print 'len gap:', len(feature)
                if "gap_type" in feature.qualifiers:
                    if feature.qualifiers["gap_type"][0] == "within scaffold":
                        continue
                    else:
                        # print 'new gap type:', feature.qualifiers["gap_type"][0]
                        pass

                if len(feature) != 200:
                    continue
                stop = feature.location.start -1
                sub_record = record[start:stop]
                start = feature.location.end+1
                if len(sub_record)>=50000: # 10000
                    big_contig_length+=len(sub_record)
                    long_record_list.append(sub_record)

                # count number of cds
                n_CDS = 0
                for feature in sub_record.features:
                    if feature.type == 'CDS':
                        n_CDS+=1
                if n_CDS == 0:
                    accession2n_contigs_without_cds[accession] += 1
                else:
                    # count number of BBH Chlamydiae
                    n_BBH_Chlamydiae = 0
                    for feature in sub_record.features:
                        if feature.type == 'CDS':
                            if not 'pseudo' in feature.qualifiers:
                                if feature.qualifiers['locus_tag'][0] in locus_BBH_chlamydiae:
                                    n_BBH_Chlamydiae+=1
                    if n_BBH_Chlamydiae == 0:
                        accession2n_countigs_without_BBH_chlamydiae[accession] += 1

        stop = len(record)
        sub_record = record[start:stop]
        if len(sub_record)>=10000:
            big_contig_length+=len(sub_record)
            long_record_list.append(sub_record)
        #print accession, 'number of big contigs:', len(long_record_list)

        len_coding = 0
        rrna_count_16 = 0
        rrna_count_23 = 0
        rrna_count_5 = 0
        trna_count = 0
        # if small plasmid < 10kb
        if len(long_record_list) == 0:
            long_record_list = [record]
            big_contig_length = len(record.seq)
        for one_record in long_record_list:
            #print one_record
            #len_seq = len(str(one_record.seq).replace("N", ""))
            for n, feature in enumerate(one_record.features):
                if feature.type in ['rRNA','tRNA','CDS'] and not 'pseudo' in feature.qualifiers:
                    len_coding += len(feature)
                if feature.type == 'tRNA':
                    trna_count+=1
                if feature.type == 'rRNA':
                    if 'product' in feature.qualifiers:
                        if '16S' in feature.qualifiers['product'][0]:
                            rrna_count_16+=1
                        if '23S' in feature.qualifiers['product'][0]:
                            rrna_count_23+=1
                        if '5S' in feature.qualifiers['product'][0]:
                            rrna_count_5+=1
                    else:
                        print("Unknown rRNA type:", feature)
        accession2density[accession] = round((len_coding/float(big_contig_length))*100,2)
        accession2n_trna[accession] = trna_count
        accession2n_rrna[accession] = [rrna_count_16, rrna_count_23, rrna_count_5]
        accession2big_contig_length[accession] = big_contig_length
        #print round((len_coding/float(big_contig_length))*100,2), trna_count, rrna_count_16, rrna_count_23, rrna_count_5
    return accession2density, accession2n_trna, accession2n_rrna, accession2big_contig_length, \
           accession2n_countigs_without_BBH_chlamydiae, accession2n_contigs_without_cds



def taxonomical_form(db_name, hit_number=1):

    '''
    contruct of an html form with accessions, superkingdom and phylum selection
    possibility to nest more levels

    '''

    server, db = manipulate_biosqldb.load_db(db_name)

    sql1 = 'select accession, bioentry.description from biosqldb.bioentry' \
           ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id' \
           ' where biodatabase.name = "%s"' % db_name

    accessions = server.adaptor.execute_and_fetchall(sql1,)

    #accessions = [accession + ('-',) for accession in accessions]


    form = '<p><label for="id_genome_taxonomy">Genome:</label><select id="genome" name="Genome">\n'
    form+='    <option value="">--</option>\n'
    for accession in accessions:
        form+='<option value="%s">%s</option>\n' % (accession[0], accession[1])

    form+='</select></p>\n'
    #print form

    # pour chaque accession, faire un choix avec les superkingdom

    form+='<p><label for="id_superkingdom_taxonomy">Superkingdom:</label><select id="superkingdom" name="Superkingdom">\n'
    form+='    <option value="">--</option>\n'
    for accession in accessions:

        sql2 = 'select t3.superkingdom ' \
              ' from blastnr_blastnr_hits_%s as t1 ' \
              ' inner join blastnr_blastnr_hits_taxonomy_filtered_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
              ' inner join blastnr_blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
              '  where t1.hit_number=%s group by t3.superkingdom' % (db_name, accession[0],db_name, accession[0], hit_number)

        #print sql2

        superkinkdoms = [i[0] for i in server.adaptor.execute_and_fetchall(sql2,)]

        for superkinkdom in superkinkdoms:
            if superkinkdom == "-":
                #print superkinkdoms, accession
                form+='    <option value="%s_%s" class="%s">%s</option>\n' % (accession[0], 'unclassified', accession[0], 'unclassified')
            else:
                form+='    <option value="%s_%s" class="%s">%s</option>\n' % (accession[0], superkinkdom, accession[0], superkinkdom)
    form+='</select></p>\n'



    form+='<p><label for="id_phylum_taxonomy">Phylum:</label><select id="phylum" name="Phylum">\n'
    form+='    <option value="">--</option>\n'
    for accession in accessions:

        sql3 = 'select distinct t3.superkingdom, t3.phylum ' \
              ' from blastnr_blastnr_hits_%s as t1 ' \
              ' inner join blastnr_blastnr_hits_taxonomy_filtered_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
              ' inner join blastnr_blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
              ' inner join blastnr_blastnr_hsps_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
              ' where t1.hit_number=%s' % (db_name, accession[0],db_name, accession[0], db_name, accession[0], hit_number)
        #print sql3
        phylums = server.adaptor.execute_and_fetchall(sql3,)

        #print sql3

        for phylum in phylums:
            #print accession, phylum
            if phylum[0] == '-' and phylum[1] != '-':
                #print 'cas 1', accession
                form+='    <option value="%s_%s_%s" class="%s_%s">%s</option>\n' % (accession[0], 'unclassified', phylum[1], accession[0],'unclassified', phylum[1])
            elif phylum[0] == '-' and phylum[1] == '-':
                #print 'cas 2', accession

                form+='    <option value="%s_%s_%s" class="%s_%s">%s</option>\n' % (accession[0], 'unclassified', 'unclassified', accession[0],'unclassified', 'unclassified')
            elif phylum[0] != '-' and phylum[1] == '-':
                #print 'cas 3', accession
                form+='    <option value="%s_%s_%s" class="%s_%s">%s</option>\n' % (accession[0], phylum[0], 'unclassified', accession[0], phylum[0], 'unclassified')



            else:
                form+='    <option value="%s_%s_%s" class="%s_%s">%s</option>\n' % (accession[0], phylum[0], phylum[1], accession[0], phylum[0], phylum[1])
    form+='</select></p>\n'
    return form

def clean_multispecies_blastnr_record(db_name, create_new_sql_tables = False):


    server, db = manipulate_biosqldb.load_db(db_name)
    sql = 'select accession from bioentry' \
      ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id' \
      ' and biodatabase.name = "%s"' % db_name

    all_accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    #all_accessions = ['Rhab']
    if create_new_sql_tables:
        for accession in all_accessions:
            sql_drop = 'DROP TABLE IF EXISTS blastnr_blastnr_hits_taxonomy_filtered_%s;' % (db_name, accession)
            sql_blast_taxonomy = ' CREATE TABLE blastnr_blastnr_hits_taxonomy_filtered_%s (nr_hit_id INT, ' \
                         ' subject_taxon_id int,' \
                         ' INDEX subject_taxon_id (subject_taxon_id))' % (db_name, accession)

            sql_blast_taxonomy2 = 'ALTER TABLE blastnr_blastnr_hits_taxonomy_filtered_%s ADD CONSTRAINT fk_blast_hit_id_filter_%s ' \
                                  'FOREIGN KEY (nr_hit_id) REFERENCES blastnr_blastnr_hits_%s(nr_hit_id);' % (db_name, accession, accession, db_name, accession)

            #print sql_drop
            server.adaptor.execute(sql_drop,)
            server.adaptor.commit()
            #print sql_blast_taxonomy
            server.adaptor.execute(sql_blast_taxonomy,)
            server.adaptor.commit()
            #print sql_blast_taxonomy2
            server.adaptor.execute(sql_blast_taxonomy2,)
            server.adaptor.commit()


    total_ambiguous_taxons = 0
    for accession in all_accessions:
        #print accession, total_ambiguous_taxons
        sql = 'select nr_hit_id from blastnr_blastnr_hits_%s as t1' % (db_name, accession) #  where hit_number=1

        all_blastnr_hit_ids = server.adaptor.execute_and_fetchall(sql,)

        more_than_one_genus_path = 0
        for blast_id in all_blastnr_hit_ids:
            blast_id = int(blast_id[0])

            sql = 'select  t3.locus_tag, t1.nr_hit_id, t1.subject_taxon_id,  t2.phylum, t2.order, t2.family, t2.genus, t2.species, t2.superkingdom, t3.hit_number, t3.subject_title from blastnr_blastnr_hits_taxonomy_%s as t1 ' \
                  ' inner join blastnr_blastnr_hits_%s as t3 on t1.nr_hit_id=t3.nr_hit_id ' \
                  ' inner join blastnr_blastnr_taxonomy as t2 on t1.subject_taxon_id=t2.taxon_id ' \
                  ' where t1.nr_hit_id=%s' % (db_name, accession, db_name, accession, blast_id)
            one_blastnr_taxons = server.adaptor.execute_and_fetchall(sql,)
            #print 'blast_id', blast_id, 'n taxons', len(one_blastnr_taxons)

            # build dictionnary
            one_nr_hit_taxons = {}
            for i in one_blastnr_taxons:
                one_nr_hit_taxons[i[2]] = {
                    'query_locus_tag': i[0],
                    'nr_hit_id': i[1],
                    'phylum': i[3],
                    'order': i[4],
                    'family': i[5],
                    'genus': i[6],
                    'species': i[7],
                    'superkingdom': i[8],
                    'hit_number': i[9],
                    'description': i[10]
                }

            '''
            best_classified_path_length = 0
            for one_taxon in one_blastnr_taxons:

                # flag to check if redundancy was not already observed in the taxon group
                temp_count = 0
                for one_rank in one_taxon[3:]:
                    if one_rank != '-':
                        temp_count += 1
                if temp_count > best_classified_path_length:
                    best_classified_path_length = temp_count

            # equally lon path extraction
            '''
            unique_taxon_paths = {}
            poor_path = {}
            # iter dictionnary
            for one_subject_taxon_id in one_nr_hit_taxons:

                # removing data from the phage (small contig identical to the phage)
                phage = ['Rhab_01766','Rhab_01767','Rhab_01768','Rhab_01769','Rhab_01770','Rhab_01771','Rhab_01772']
                #print one_nr_hit_taxons[one_subject_taxon_id]
                if one_nr_hit_taxons[one_subject_taxon_id]['query_locus_tag'] in phage:
                    break

                # check if phylum, family, genus, species are not null
                temp_count = 0
                for one_rank in ['phylum', 'order', 'family', 'genus', 'species', 'superkingdom']:
                    if one_nr_hit_taxons[one_subject_taxon_id][one_rank] != '-':
                        temp_count += 1
                one_nr_hit_taxons[one_subject_taxon_id]['path_completeness'] = temp_count


                # check if an idenical subject taxon path is already present in equall_taxon_path
                # different taxon id can have identical phylum-genus path: i.e strains of the same species, subspecies,..
                if one_subject_taxon_id not in unique_taxon_paths:

                    counted = False
                    # check if taxon path with same family and same genus is present
                    # check if taxon path with same family and no genus ('-') is present
                    ### => replace it if genus exist in the new path
                    for unique_taxon_path in unique_taxon_paths.keys():
                        if (unique_taxon_paths[unique_taxon_path]['family'] == one_nr_hit_taxons[one_subject_taxon_id]['family']
                            and unique_taxon_paths[unique_taxon_path]['genus'] == one_nr_hit_taxons[one_subject_taxon_id]['genus']) \
                                or (unique_taxon_paths[unique_taxon_path]['family'] == one_nr_hit_taxons[one_subject_taxon_id]['family']
                                    and one_nr_hit_taxons[one_subject_taxon_id]['genus'] == '-'):
                            unique_taxon_paths[unique_taxon_path]['taxon_path_frequency'] += 1
                            counted = True
                            # one match, break the loop
                            break
                        # if same family as the reference but the reference has no genus and the new taxon has one
                        elif (unique_taxon_paths[unique_taxon_path]['family'] == one_nr_hit_taxons[one_subject_taxon_id]['family']
                                and unique_taxon_paths[unique_taxon_path]['genus'] == '-'):

                            # retrieve the count of the uncomplete path, replace it and delete it
                            unique_taxon_paths[one_subject_taxon_id] = one_nr_hit_taxons[one_subject_taxon_id]
                            unique_taxon_paths[one_subject_taxon_id]['taxon_path_frequency'] = unique_taxon_paths[unique_taxon_path]['taxon_path_frequency'] + 1
                            del unique_taxon_paths[unique_taxon_path]
                            counted = True
                            # one match break the loop
                            break
                    # if no match with any unique path, create a new unique entry with count 1 (new family / genus combination)
                    if not counted:
                        unique_taxon_paths[one_subject_taxon_id] = one_nr_hit_taxons[one_subject_taxon_id]
                        unique_taxon_paths[one_subject_taxon_id]['taxon_path_frequency'] = 1
                # if two time identical taxon (can not happen, oder?)
                else:
                    unique_taxon_paths[one_subject_taxon_id]['taxon_path_frequency']  += 1

            '''
            # count the most frequent paths and keep all equally frequent path
            if len(equall_taxon_path) > 1:
                max_frequency = 0
                for i in equall_taxon_path:
                    if equall_taxon_path[i][0] > max_frequency:
                        max_frequency = equall_taxon_path[i][0]

                keep_frequency = []
                for i in equall_taxon_path:
                    if equall_taxon_path[i][0] == max_frequency:
                        keep_frequency.append(equall_taxon_path[i])


                for i in keep_frequency:
                    print 'frequency:', i[0], 'path length', i[2], i[1][2:]
            '''



            poor_path = {}
            ok_path = {}
            # if with have path with different family/genus
            # check if all some are not complete (if only 2 rank/5 are recorded) ('chlamydiae', 'chlamydiales', '-', '-', '-', '-')
            if len(unique_taxon_paths) > 1:

                for one_path_taxid in unique_taxon_paths:
                    if unique_taxon_paths[one_path_taxid]['path_completeness'] < 4:
                        #print 'poor path!', unique_taxon_paths[one_path_taxid]
                        poor_path[one_path_taxid] = unique_taxon_paths[one_path_taxid]
                    else:
                        ok_path[one_path_taxid] = unique_taxon_paths[one_path_taxid]
                # if we have at least one good path

                if len(ok_path) > 1:

                    # hit associated with different families/genus
                    more_than_one_genus_path += 1
                    for taxon_id in ok_path:

                        sql = 'INSERT INTO blastnr_blastnr_hits_taxonomy_filtered_%s (' \
                        'nr_hit_id, subject_taxon_id) values (%s, %s)'  % (db_name,
                                                                           accession,
                                                                           ok_path[taxon_id]['nr_hit_id'],
                                                                           taxon_id)

                        server.adaptor.execute(sql,)
                        server.adaptor.commit()
                        '''
                        print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (len(ok_path),
                                    ok_path[taxon_id]['nr_hit_id'],
                                    ok_path[taxon_id]['hit_number'],
                                    ok_path[taxon_id]['description'],
                                    ok_path[taxon_id]['query_locus_tag'],
                                    taxon_id,
                                    ok_path[taxon_id]['superkingdom'],
                                    ok_path[taxon_id]['phylum'],
                                    ok_path[taxon_id]['order'],
                                    ok_path[taxon_id]['family'],
                                    ok_path[taxon_id]['genus'])
                        '''

                # if we do not have any good path
                elif len(ok_path) == 0 and len(poor_path) == 0:
                    #print 'No path: Impossible!'
                    import sys
                    sys.exit()
                elif len(ok_path) == 1:

                    for taxon_id in ok_path:

                        sql = 'INSERT INTO blastnr_blastnr_hits_taxonomy_filtered_%s (' \
                        'nr_hit_id, subject_taxon_id) values (%s, %s)'  % (db_name,
                                                                           accession,
                                                                           ok_path[taxon_id]['nr_hit_id'],
                                                                           taxon_id)

                        server.adaptor.execute(sql,)
                        server.adaptor.commit()
                else:
                    # no good path
                    for taxon_id in poor_path:

                        sql = 'INSERT INTO blastnr_blastnr_hits_taxonomy_filtered_%s (' \
                        'nr_hit_id, subject_taxon_id) values (%s, %s)'  % (db_name,
                                                                           accession,
                                                                           poor_path[taxon_id]['nr_hit_id'],
                                                                           taxon_id)

                        server.adaptor.execute(sql,)
                        server.adaptor.commit()

            # only one path
            else:
                for taxon_id in unique_taxon_paths:

                    sql = 'INSERT INTO blastnr_blastnr_hits_taxonomy_filtered_%s (' \
                        'nr_hit_id, subject_taxon_id) values (%s, %s)'  % (db_name,
                                                                           accession,
                                                                           unique_taxon_paths[taxon_id]['nr_hit_id'],
                                                                           taxon_id)

                    server.adaptor.execute(sql,)
                    server.adaptor.commit()


    #print "total_ambiguous_taxons", total_ambiguous_taxons

    '''
    nr_hit_id |
    subject_taxon_id |
    taxon_id |

    no_rank            |
    superkingdom |
    kingdom |
    subkingdom |
    superphylum |
    phylum         |
    subphylum |
    superclass |
    class               |
    subclass |
    superorder |
    order           |
    suborder |
    superfamily |
    family           |
    subfamily |
    genus         |
    subgenus |
    species |
    species_subgroup |
    species_group |
    subspecies |
    tribe |
    infraorder |
    subtribe |
    forma |
    infraclass |
    varietas |
    parvorder |
    '''



def get_best_hit_excluding_one_family():
    sql='select t3.subject_taxon_id,t4.no_rank, t4.phylum, t4.order, t4.family, B.* from (select * from biosqldb.orthology_detail_chlamydia_03_15 as t1 where t1.accession="NC_015713") A left join (select t2.nr_hit_id,t2.locus_tag,t2.subject_kingdom,t2.subject_accession from blastnr.blastnr_hits_chlamydia_03_15_NC_015713 as t2) B on A.locus_tag=B.locus_tag inner join blastnr.blastnr_hits_taxonomy_chlamydia_03_15_NC_015713 as t3 on B.nr_hit_id=t3.nr_hit_id inner join blastnr_blastnr_taxonomy as t4 on t3.subject_taxon_id=t4.taxon_id where t4.family !="Simkaniaceae";'

def locus_tag2n_nr_hits(db_name, genome_accession, exclude_family = False):
    #print 'exclude family', exclude_family
    server, db = manipulate_biosqldb.load_db(db_name)
    if not exclude_family:
        sql = 'select  A.locus_tag, B.n_hits from (select * from orthology_detail as t1 where t1.accession="%s") A ' \
              ' left join (select t2.nr_hit_id,t2.locus_tag,t2.subject_kingdom,t2.subject_accession,count(t2.locus_tag) as n_hits ' \
              ' from blastnr_blastnr_hits_%s as t2 ' \
              ' group by locus_tag) B on A.locus_tag=B.locus_tag;' %(db_name, genome_accession, db_name, genome_accession)

        sql = 'select locus_tag, count(*) from custom_tables_locus2seqfeature_id t1 ' \
              ' left join blastnr_blastnr as t2 on t1.seqfeature_id=t2.seqfeature_id ' \
              ' inner join biosqldb.bioentry as t3 on t2.query_bioentry_id=t3.bioentry_id ' \
              ' where t3.accession="%s" group by locus_tag;' % (db_name, db_name,genome_accession)
    else:
        sql = 'select  A.locus_tag, B.n_hits from (select * from orthology_detail as t1 where t1.accession="%s") A ' \
              ' left join (select t2.nr_hit_id,t2.locus_tag,t2.subject_kingdom,t2.subject_accession, count(t2.locus_tag) as n_hits ' \
              ' from blastnr_blastnr_hits_%s as t2 inner join blastnr_blastnr_hits_taxonomy_filtered_%s as t3 on t3.nr_hit_id=t2.nr_hit_id' \
              ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id where t4.family!="%s"' \
              ' group by locus_tag) B on A.locus_tag=B.locus_tag;' % (db_name,
                                                                     genome_accession,
                                                                     db_name,
                                                                     genome_accession,
                                                                     db_name,
                                                                     genome_accession,
                                                                     exclude_family)


    #print sql
    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

def collect_genome_statistics(biodb, sqlite=False):
    from Bio.SeqUtils import GC
    server, db = manipulate_biosqldb.load_db(biodb,sqlite=sqlite)

    sql = 'select accession, seq as length from bioentry as t1 ' \
          ' inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id ' \
          ' inner join biosequence as t3 on t1.bioentry_id=t3.bioentry_id where t2.name="%s";' % biodb
    #print 'getting accession2genome_sequence...'

    accession2genome_sequence = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select t3.accession,t3.description, count(*) from annotation_seqfeature_id2locus t1 ' \
           ' inner join annotation_seqfeature_id2CDS_annotation t2 ' \
           ' on t1.seqfeature_id=t2.seqfeature_id inner join bioentry t3 on t1.bioentry_id=t3.bioentry_id ' \
           ' inner join biodatabase t4 on t3.biodatabase_id=t4.biodatabase_id ' \
           ' where t4.name="%s" group by t1.bioentry_id;' % (biodb, biodb, biodb)


    #print 'getting accession2genome_data...'
    # organism name, protein encoding ORF number
    #print sql
    accession2genome_data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    accession2density, \
    accession2n_trna, \
    accession2n_rrna, \
    accession2big_contig_length, \
    accession2n_countigs_without_BBH_chlamydiae, \
    accession2n_contigs_without_cds = accession2coding_density(biodb, sqlite=sqlite)

    #print '<td colspan="6"><table width="800" border=0  class=table_genomes>'

    #print 'preparing data...'
    genomes_data = []

    for accession in accession2genome_data:
        accession = accession.split('.')[0]
        #print accession
        one_genome_data = []
        one_genome_data.append(accession)
        try:
            one_genome_data.append(round(GC(accession2genome_sequence[accession].replace("N", "")),2))
        except:
            one_genome_data.append('-')
        one_genome_data.append(accession2genome_data[accession][1])
        try:
            one_genome_data.append(accession2genome_sequence[accession].count(200*'N') + 1)
        except:
            one_genome_data.append('-')
        try:
            size = len(accession2genome_sequence[accession]) - (200*accession2genome_sequence[accession].count(200*'N'))
            one_genome_data.append(str(size))
        except:
            one_genome_data.append('-')
        one_genome_data.append(accession2genome_data[accession][0])
        genomes_data.append(one_genome_data)

    #['NC_017080', 73.29, 2955L, 1, '3803225', 'Phycisphaera mikurensis NBRC 102666'], ['LVAZ01000000', 67.26, 3793L, 662, '4951698', 'Victivallales bacterium Lenti_02']
    # accession, GC, n_CDS, n_contigs, genome_size, description
    sql = 'create table if not EXISTS genomes_info_%s (ACCESSION VARCHAR(200), ' \
          ' GC FLOAT, n_CDS INT ,n_contigs INT, genome_size INT, ' \
          ' n_tRNA INT, n_16S INT, n_23S INT, n_5S INT, percent_non_coding FLOAT, big_contig_length INT, ' \
          ' n_no_CDS INT,' \
          ' n_no_BBH_chlamydiae INT,' \
          ' description VARCHAR (20000))' % biodb

    server.adaptor.execute_and_fetchall(sql,)

    for one_genome in genomes_data:
        sql = 'insert into genomes_info_%s values ("%s", %s, %s, %s, %s, %s, %s, %s, %s, %s,%s,%s,%s,"%s")' % (biodb,
                                                                           one_genome[0],
                                                                           one_genome[1],
                                                                           one_genome[2],
                                                                           one_genome[3],
                                                                           one_genome[4],
                                                                           accession2n_trna[one_genome[0]],
                                                                           accession2n_rrna[one_genome[0]][0],
                                                                           accession2n_rrna[one_genome[0]][1],
                                                                           accession2n_rrna[one_genome[0]][2],
                                                                           100-accession2density[one_genome[0]],
                                                                           accession2big_contig_length[one_genome[0]],
                                                                           accession2n_contigs_without_cds[one_genome[0]],
                                                                           accession2n_countigs_without_BBH_chlamydiae[one_genome[0]],
                                                                           one_genome[5])
        #print sql
        server.adaptor.execute(sql,)
        server.commit()

def get_comparative_subtable(biodb,
                             table_name,
                             first_col_name,
                             taxon_list,
                             exclude_taxon_list,
                             ratio=1,
                             single_copy=False,
                             accessions=False,
                             cache=False):

    import pandas
    import numpy

    '''
    ratio: ratio of missing data tolaration
    single_copy: only consider singly copy locus
    accessions (bol): compare accession and not taxons (differentiate plasmids from chromosome)
    '''
    
    import MySQLdb
    import os
    import pandas
    mysql_host = 'localhost'
    mysql_user = 'root'
    mysql_pwd = os.environ['SQLPSW']
    mysql_db = 'comparative_tables'
    conn = MySQLdb.connect(host=mysql_host,
                                user=mysql_user,
                                passwd=mysql_pwd,
                                db=mysql_db)
    
    if not accessions:
        cache_id = f"{biodb}_{table_name}2"   
        count_df = cache.get(cache_id)

        if not isinstance(count_df, pandas.DataFrame):
            sql = f'select * from comparative_tables.{table_name}_{biodb}'
            count_df = pandas.read_sql(sql, conn, index_col=first_col_name)
            cache.set(f"{biodb}_{table_name}", count_df)
    else:
        count_df = cache.get(f"{biodb}_{table_name}_acc")

        if not count_df:
            sql = f'select * from comparative_tables.{table_name}_accessions_{biodb}'
            count_df = pandas.read_sql(sql, conn,index_col=first_col_name)
            cache.set(f"{biodb}_{table_name}_acc", count_df)      
    
     # convert to integer
    count_df = count_df.apply(pandas.to_numeric, args=('coerce',))
    # rename columns (interger are not ok for df.query())
    count_df.columns = ["taxid_%s" % i for i in list(count_df)]

    include_cols = ["taxid_%s" % i for i in taxon_list]
    if len(exclude_taxon_list) > 0:
        # apply exclude filter
        exclude_string = ' & '.join(["taxid_%s==0" % i for i in exclude_taxon_list])       
        count_df2 = count_df.query(exclude_string).loc[:,include_cols]
    else:
        count_df2 = count_df.loc[:,include_cols]
  
    # calculate limit when accepting missing data
    limit = len(taxon_list)*ratio
    
    if not single_copy: 
        
        return count_df2[(count_df2 > 0).sum(axis=1) >= limit], count_df
    else:
        # identify and remove rows with paralogs
        groups_with_paralogs = count_df2[(count_df2 > 1).sum(axis=1) > 0].index
        count_df2 = count_df2.drop(groups_with_paralogs)

        return count_df2[(count_df2 == 1).sum(axis=1) >= limit], count_df       

    
def best_hit_classification(db_name, accession):
    sql = 'select A.*, t4.superkingdom, t4.kingdom, t4.phylum, t4.order, t4.family, t4.genus, t4.species' \
          ' from (select locus_tag, TM, SP, gene, product from orthology_detail as t1 where t1.accession="%s" group by locus_tag) A' \
          ' left join (select t2.nr_hit_id, locus_tag from blastnr_blastnr_hits_%s as t2  where hit_number=1) B on A.locus_tag=B.locus_tag' \
          ' left join blastnr_blastnr_hits_taxonomy_filtered_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
          ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id;' % (db_name,
                                                                                                                        accession,
                                                                                                                        db_name,
                                                                                                                        accession,                                                                                                              db_name,
                                                                                                                        accession)

    import pandas
    #print sql
    server, db = manipulate_biosqldb.load_db(db_name)

    data = [i for i in server.adaptor.execute_and_fetchall(sql)]
    table = pandas.DataFrame(data, columns=['locus_tag', 'TM', 'SP', 'gene', 'product', 'superkingdom','kingdom','phylum','order','family','genus','species'])
    return table

def best_hit_phylum_and_protein_length(db_name, accession):
    sql = 'select A.locus_tag, char_length(A.translation), t4.superkingdom' \
          ' from (select locus_tag, translation from orthology_detail as t1 where t1.accession="%s" group by locus_tag) A' \
          ' left join (select t2.nr_hit_id, locus_tag from blastnr_blastnr_hits_%s as t2  where hit_number=1) B on A.locus_tag=B.locus_tag' \
          ' left join blastnr_blastnr_hits_taxonomy_filtered_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
          ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id;' % (db_name,
                                                                                                                        accession,
                                                                                                                        db_name,
                                                                                                                        accession,                                                                                                              db_name,
                                                                                                                        accession)
    #print sql
    import pandas
    server, db = manipulate_biosqldb.load_db(db_name)

    data = [i for i in server.adaptor.execute_and_fetchall(sql)]
    table = pandas.DataFrame(data, columns=['locus_tag', 'protein_length', 'superkingdom'])
    return table



def locus_tag2n_blast_superkingdom(db_name, accession, superkingdom="Bacteria", exclude_family=False):
    if not exclude_family:
        sql = 'select A.locus_tag, count(*)' \
              ' from (select locus_tag, translation from orthology_detail as t1 where t1.accession="%s" group by locus_tag) A' \
              ' left join (select t2.nr_hit_id, locus_tag from blastnr_blastnr_hits_%s as t2) B on A.locus_tag=B.locus_tag' \
              ' left join blastnr_blastnr_hits_taxonomy_filtered_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
              ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
              ' where t4.superkingdom="%s" group by locus_tag;' % (db_name,
                                                                        accession,
                                                                        db_name,
                                                                        accession,                                                                                                              db_name,
                                                                        accession,
                                                                        superkingdom)
        sql = 'select locus_tag, count(*) from custom_tables_locus2seqfeature_id t1 ' \
              ' inner join blastnr_blastnr as t2 on t1.seqfeature_id=t2.seqfeature_id ' \
              ' inner join biosqldb.bioentry as t3 on t2.query_bioentry_id=t3.bioentry_id ' \
              ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t2.subject_taxid ' \
              'where t3.accession="%s" and t4.superkingdom="%s" group by locus_tag;' % (db_name,
                                                                                        db_name,
                                                                                        accession,
                                                                                        superkingdom)
    else:

        sql = 'select A.locus_tag, count(*)' \
              ' from (select locus_tag, translation from orthology_detail as t1 where t1.accession="%s" group by locus_tag) A' \
              ' left join (select t2.nr_hit_id, locus_tag from blastnr_blastnr_hits_%s as t2) B on A.locus_tag=B.locus_tag' \
              ' left join blastnr_blastnr_hits_taxonomy_filtered_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
              ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
              ' where t4.superkingdom="%s" and t4.family!="%s" group by locus_tag;' % (db_name,
                                                                        accession,
                                                                        db_name,
                                                                        accession,                                                                                                              db_name,
                                                                        accession,
                                                                        superkingdom,
                                                                        exclude_family)


    #print sql
    server, db = manipulate_biosqldb.load_db(db_name)
    data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    return data

def locus_tag2n_blast_bacterial_phylum(db_name, accession, phylum="Chlamydiae", reverse=False, exclude_family = False):
    if not reverse:
        if not exclude_family:
            sql = 'select A.locus_tag, count(*)' \
                  ' from (select locus_tag, translation from orthology_detail as t1 where t1.accession="%s" group by locus_tag) A' \
                  ' left join (select t2.nr_hit_id, locus_tag from blastnr_blastnr_hits_%s as t2) B on A.locus_tag=B.locus_tag' \
                  ' left join blastnr_blastnr_hits_taxonomy_filtered_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
                  ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
                  ' where t4.superkingdom="Bacteria" and t4.phylum="%s" group by locus_tag;' % (db_name,
                                                                            accession,
                                                                            db_name,
                                                                            accession,                                                                                                              db_name,
                                                                            accession,
                                                                            phylum)
            sql = 'select locus_tag, count(*) from custom_tables_locus2seqfeature_id t1 ' \
                  ' inner join blastnr_blastnr as t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                  ' inner join biosqldb.bioentry as t3 on t2.query_bioentry_id=t3.bioentry_id ' \
                  ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t2.subject_taxid ' \
                  'where t3.accession="%s" and t4.phylum="%s" group by locus_tag;' % (db_name,
                                                                                        db_name,
                                                                                        accession,
                                                                                        phylum)
        else:
            sql = 'select A.locus_tag, count(*)' \
                  ' from (select locus_tag, translation from orthology_detail as t1 where t1.accession="%s" group by locus_tag) A' \
                  ' left join (select t2.nr_hit_id, locus_tag from blastnr_blastnr_hits_%s as t2) B on A.locus_tag=B.locus_tag' \
                  ' left join blastnr_blastnr_hits_taxonomy_filtered_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
                  ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
                  ' where t4.superkingdom="Bacteria" and t4.phylum="%s" and t4.family!="%s" group by locus_tag;' % (db_name,
                                                                            accession,
                                                                            db_name,
                                                                            accession,                                                                                                              db_name,
                                                                            accession,
                                                                            phylum,
                                                                            exclude_family)
    else:
        if not exclude_family:
            sql = 'select A.locus_tag, count(*)' \
                  ' from (select locus_tag, translation from orthology_detail as t1 where t1.accession="%s" group by locus_tag) A' \
                  ' left join (select t2.nr_hit_id, locus_tag from blastnr_blastnr_hits_%s as t2) B on A.locus_tag=B.locus_tag' \
                  ' left join blastnr_blastnr_hits_taxonomy_filtered_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
                  ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
                  ' where t4.superkingdom="Bacteria" and t4.phylum !="%s" group by locus_tag;' % (db_name,
                                                                            accession,
                                                                            db_name,
                                                                            accession,                                                                                                              db_name,
                                                                            accession,
                                                                            phylum)

            sql = 'select locus_tag, count(*) from custom_tables_locus2seqfeature_id t1 ' \
                  ' inner join blastnr_blastnr as t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                  ' inner join biosqldb.bioentry as t3 on t2.query_bioentry_id=t3.bioentry_id ' \
                  ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t2.subject_taxid ' \
                  'where t3.accession="%s" and t4.superkingdom="Bacteria" and t4.phylum !="%s" group by locus_tag;' % (db_name,
                                                                                        db_name,
                                                                                        accession,
                                                                                        phylum)
        else:
            sql = 'select A.locus_tag, count(*)' \
                  ' from (select locus_tag, translation from orthology_detail as t1 where t1.accession="%s" group by locus_tag) A' \
                  ' left join (select t2.nr_hit_id, locus_tag from blastnr_blastnr_hits_%s as t2) B on A.locus_tag=B.locus_tag' \
                  ' left join blastnr_blastnr_hits_taxonomy_filtered_%s as t3 on t3.nr_hit_id=B.nr_hit_id' \
                  ' left join blastnr_blastnr_taxonomy as t4 on t4.taxon_id = t3.subject_taxon_id' \
                  ' where t4.superkingdom="Bacteria" and t4.phylum !="%s" and t4.family!="%s" group by locus_tag;' % (db_name,
                                                                            accession,
                                                                            db_name,
                                                                            accession,                                                                                                              db_name,
                                                                            accession,
                                                                            phylum, exclude_family)

    #print sql
    server, db = manipulate_biosqldb.load_db(db_name)
    data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    return data







def locus_tag2orthogroup(db_name):
    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select locus_tag, orthogroup from orthology_detail group by locus_tag, orthogroup' % db_name

    return manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


def locus_tag2presence_in_n_genomes(db_name):
    '''
    return a dictionnary of all locus tag and the number of genomes in which they have one or multiple homolog(s)
    '''

    from chlamdb.biosqldb import mysqldb_load_mcl_output
    server, db = manipulate_biosqldb.load_db(db_name)
    orthogroup2family_size = mysqldb_load_mcl_output.get_family_size(server, db_name)
    locus_tag2orthogroup_dico = locus_tag2orthogroup(db_name)

    locus_tag2presence_absence = {}
    for locus in locus_tag2orthogroup_dico:
        locus_tag2presence_absence[locus] = orthogroup2family_size[locus_tag2orthogroup_dico[locus]]
    return locus_tag2presence_absence

def locus_tag2n_paralogs(db_name, genome_accession):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select orthogroup, count(*) as paralogs ' \
          ' from orthology_detail ' \
          ' where accession="%s" group by orthogroup order by paralogs DESC;' % (db_name, genome_accession)

    return server.adaptor.execute_and_fetchall(sql,)

def orthogroup2gene(db_name, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    if not accession:
        sql =  'select orthogroup, gene from orthology_detail' % db_name
    else:
        sql = 'select orthogroup, gene from orthology_detail where accession="%s"' % (db_name, accession)


    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2gene = {}

    for row in data:
        group = row[0]
        gene = row[1].split("_")[0]
        if group not in ortho2gene:
            ortho2gene[group] = {}
            ortho2gene[group][gene] = 1
        else:
            if gene in ortho2gene[group]:
                ortho2gene[group][gene] += 1
            else:
                ortho2gene[group][gene] = 1

    return ortho2gene


def orthogroup2cog(db_name, accession=False, group_list=False): # group_list,

    server, db = manipulate_biosqldb.load_db(db_name)

    #group_list_form = '"' + '","'.join(group_list) + '"'
    '''
    if not accession:
        sql =  'select t1.orthogroup, t2.cog_id from (select * from orthology_detail ' \
               ' where orthogroup in (%s)) t1 left join COG_locus_tag2gi_hit ' \
               'as t2 on t1.locus_tag=t2.locus_tag;' % (db_name, group_list_form, db_name)
    else:
        sql = 'select t1.orthogroup, t2.cog_id from orthology_detail as t1 left join COG_locus_tag2gi_hit ' \
              'as t2 on t1.locus_tag=t2.locus_tag' % (db_name, db_name)
    '''
    if not accession:
        sql =  'select orthogroup_name,COG_name from COG_seqfeature_id2best_COG_hit t1 inner join COG_cog_names_2014 t2 on t1.hit_COG_id=t2.COG_id inner join orthology.seqfeature_id2orthogroup_%s t3 on t1.seqfeature_id=t3.seqfeature_id inner join orthology.orthogroup_%s t4 on t3.orthogroup_id=t4.orthogroup_id' % (db_name, db_name, db_name)
        sql = 'select orthogroup_name,COG_name from COG_seqfeature_id2best_COG_hit t1 inner join COG_cog_names_2014 t2 on t1.hit_COG_id=t2.COG_id inner join orthology.seqfeature_id2orthogroup_%s t3 on t1.seqfeature_id=t3.seqfeature_id inner join orthology.orthogroup_%s t4 on t3.orthogroup_id=t4.orthogroup_id' % (db_name, db_name, db_name)
        sql = 'select t1.orthogroup, t2.cog_id from orthology_detail as t1 left join COG_locus_tag2gi_hit ' \
              'as t2 on t1.locus_tag=t2.locus_tag' % (db_name, db_name)
    else:
        sql = 'select t1.orthogroup, t2.cog_id from orthology_detail as t1 left join COG_locus_tag2gi_hit ' \
              'as t2 on t1.locus_tag=t2.locus_tag' % (db_name, db_name)




    #print 'df', sql
    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2cog = {}

    for row in data:
        if row[0] not in ortho2cog:
            ortho2cog[row[0]] = {}
            ortho2cog[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2cog[row[0]]:
                ortho2cog[row[0]][row[1]] += 1
            else:
                ortho2cog[row[0]][row[1]] = 1

    return ortho2cog


def orthogroup2cog_id(db_name, accession=False): # group_list,

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select t4.orthogroup_name,hit_cog_id from COG_seqfeature_id2best_COG_hit t1 inner join COG_cog_names_2014 t2 on t1.hit_COG_id=t2.COG_id inner join orthology.seqfeature_id2orthogroup_%s t3 on t1.seqfeature_id=t3.seqfeature_id inner join orthology.orthogroup_%s t4 on t3.orthogroup_id=t4.orthogroup_id' % (db_name, db_name, db_name)

    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2cog = {}

    for row in data:
        if row[0] not in ortho2cog:
            ortho2cog[row[0]] = {}
            ortho2cog[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2cog[row[0]]:
                ortho2cog[row[0]][row[1]] += 1
            else:
                ortho2cog[row[0]][row[1]] = 1

    return ortho2cog


def orthogroup2pfam(db_name, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    if not accession:
        sql =  'select orthogroup, signature_accession from interpro ' \
               ' where analysis="Pfam"' % (db_name, )
    else:
        sql = 'select t1.orthogroup, t2.cog_id from orthology_detail as t1 left join (COG_locus_tag2gi_hit ' \
              'as t2 on t1.locus_tag=t2.locus_tag' % (db_name, db_name)

    #print 'df', sql
    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2pfam = {}

    for row in data:
        if row[0] not in ortho2pfam:
            ortho2pfam[row[0]] = {}
            ortho2pfam[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2pfam[row[0]]:
                ortho2pfam[row[0]][row[1]] += 1
            else:
                ortho2pfam[row[0]][row[1]] = 1

    return ortho2pfam


def orthogroup2pfam_id(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select A.orthogroup_name,A.signature_id from (select distinct t1.seqfeature_id,t3.orthogroup_name,t4.signature_id from interpro_interpro t1 inner join orthology.seqfeature_id2orthogroup_%s t2 on t1.seqfeature_id=t2.seqfeature_id inner join orthology.orthogroup_%s t3 on t2.orthogroup_id=t3.orthogroup_id inner join interpro_signature t4 on t1.signature_id=t4.signature_id inner join interpro_analysis t5 on t4.analysis_id=t5.analysis_id where t5.analysis_name="Pfam") A;' % (db_name, db_name,db_name)

    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2pfam = {}

    for row in data:
        if row[0] not in ortho2pfam:
            ortho2pfam[row[0]] = {}
            ortho2pfam[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2pfam[row[0]]:
                ortho2pfam[row[0]][row[1]] += 1
            else:
                ortho2pfam[row[0]][row[1]] = 1

    return ortho2pfam


def orthogroup2ko_id(db_name, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql =  'select t3.orthogroup_name,ko_id from enzyme_seqfeature_id2ko t1 inner join orthology.seqfeature_id2orthogroup_%s t2 on t1.seqfeature_id=t2.seqfeature_id inner join orthology.orthogroup_%s t3 on t2.orthogroup_id=t3.orthogroup_id;' % (db_name, db_name, db_name)

    #print 'df', sql
    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2ko = {}

    for row in data:
        if row[0] not in ortho2ko:
            ortho2ko[row[0]] = {}
            ortho2ko[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2ko[row[0]]:
                ortho2ko[row[0]][row[1]] += 1
            else:
                ortho2ko[row[0]][row[1]] = 1

    return ortho2ko


def orthogroup2interpro_id(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select A.orthogroup_name,A.interpro_id from (select distinct t1.seqfeature_id,t3.orthogroup_name,t4.interpro_id from interpro_interpro t1 inner join orthology.seqfeature_id2orthogroup_%s t2 on t1.seqfeature_id=t2.seqfeature_id inner join orthology.orthogroup_%s t3 on t2.orthogroup_id=t3.orthogroup_id inner join interpro_signature t4 on t1.signature_id=t4.signature_id) A;' % (db_name, db_name, db_name)

    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2interpro = {}

    for row in data:
        if row[0] not in ortho2interpro:
            ortho2interpro[row[0]] = {}
            ortho2interpro[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2interpro[row[0]]:
                ortho2interpro[row[0]][row[1]] += 1
            else:
                ortho2interpro[row[0]][row[1]] = 1

    return ortho2interpro


def pfam2description(db_name, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)


    if not accession:
        sql =  'select signature_accession, signature_description from interpro ' \
               ' where analysis="Pfam" group by signature_accession' % (db_name)
    else:
        sql = 'select t1.orthogroup, t2.cog_id from orthology_detail as t1 left join COG_locus_tag2gi_hit ' \
              'as t2 on t1.locus_tag=t2.locus_tag' % (db_name, db_name)


    data = server.adaptor.execute_and_fetchall(sql,)



    return manipulate_biosqldb.to_dict(data)









def orthogroup2product(db_name, accession=False):

    server, db = manipulate_biosqldb.load_db(db_name)

    if not accession:
        sql =  'select orthogroup, product from orthology_detail' % db_name
    else:
        sql = 'select orthogroup, product from orthology_detail where accession="%s"' % (db_name, accession)


    data = server.adaptor.execute_and_fetchall(sql,)

    ortho2product = {}

    for row in data:
        if row[0] not in ortho2product:
            ortho2product[row[0]] = {}
            ortho2product[row[0]][row[1]] = 1
        else:
            if row[1] in ortho2product[row[0]]:
                ortho2product[row[0]][row[1]] += 1
            else:
                ortho2product[row[0]][row[1]] = 1

    return ortho2product


'''
print 'tata'
locus_tag2best_hit_dico = locus_tag2best_hit("chlamydia_03_15", "Rhab", hit_number=1, rank=False, taxon_name=False)
print
print locus_tag2best_hit_dico.keys()[0]
print locus_tag2best_hit_dico[locus_tag2best_hit_dico.keys()[0]]

locus_tag2best_hit_euk = locus_tag2best_hit("chlamydia_03_15", "Rhab", hit_number=1, rank="superkingdom", taxon_name="Eukaryota")
print locus_tag2best_hit_euk
'''

def calculate_average_protein_identity(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select taxon_id from bioentry' \
          ' inner join biodatabase on bioentry.biodatabase_id=biodatabase.biodatabase_id where biodatabase.name="%s" group by taxon_id' % db_name

    all_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


    average_id = {}
    for i_taxon in range(0,len(all_taxons)):
        #print 'i taxon', i_taxon
        for y_taxon in range(i_taxon+1, len(all_taxons)):
            shared_groups_sql = 'select orthogroup from comparative_tables_orthology where `%s` =1 and `%s`=1' % (db_name, all_taxons[i_taxon], all_taxons[y_taxon])
            all_groups = [i[0] for i in server.adaptor.execute_and_fetchall(shared_groups_sql,)]
            identity_values = []
            #print 'y taxon', y_taxon
            for group in all_groups:
                sql_ref = 'select locus_tag from orthology_detail where taxon_id=%s and orthogroup="%s"' % (db_name, all_taxons[i_taxon], group)
                sql_query = 'select locus_tag from orthology_detail where taxon_id=%s and orthogroup="%s"' % (db_name, all_taxons[y_taxon], group)
                locus_ref = server.adaptor.execute_and_fetchall(sql_ref,)[0][0]
                locus_query = server.adaptor.execute_and_fetchall(sql_query,)[0][0]

                sql_id = 'select `%s` from orth_%s.%s where locus_tag="%s"' % (locus_ref, db_name, group, locus_query)

                id_value = server.adaptor.execute_and_fetchall(sql_id,)[0][0]

                identity_values.append(id_value)
            one_average_id = sum(identity_values) / float(len(identity_values))
            if not all_taxons[i_taxon] in average_id:
                average_id[all_taxons[i_taxon]] = {}
            average_id[all_taxons[i_taxon]][all_taxons[y_taxon]] = [one_average_id, len(identity_values)]
            #print "one_average_id",all_taxons[i_taxon],all_taxons[y_taxon], one_average_id, len(identity_values)
    #import json
    #with open('genome_identity_dico.json', 'wb') as fp:
    #    json.dump(average_id, fp)
    return average_id

def calculate_average_protein_identity_new_tables(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select taxon_id from bioentry' \
          ' inner join biodatabase on bioentry.biodatabase_id=biodatabase.biodatabase_id where biodatabase.name="%s" group by taxon_id' % db_name

    all_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


    average_id = {}
    for i_taxon in range(0,len(all_taxons)):
        #print 'i taxon', i_taxon
        for y_taxon in range(i_taxon+1, len(all_taxons)):
            shared_groups_sql = 'select orthogroup from comparative_tables_orthology where `%s` =1 and `%s`=1' % (db_name, all_taxons[i_taxon], all_taxons[y_taxon])
            all_groups = [i[0] for i in server.adaptor.execute_and_fetchall(shared_groups_sql,)]
            identity_values = []
            #print 'y taxon', y_taxon
            for i, group in enumerate(all_groups):
                if i%500 == 0:
                    print ("%s/%s" % (i, len(all_groups)))
                sql_ref = 'select locus_tag from orthology_detail where taxon_id=%s and orthogroup="%s"' % (db_name, all_taxons[i_taxon], group)
                sql_query = 'select locus_tag from orthology_detail where taxon_id=%s and orthogroup="%s"' % (db_name, all_taxons[y_taxon], group)
                locus_ref = server.adaptor.execute_and_fetchall(sql_ref,)[0][0]
                locus_query = server.adaptor.execute_and_fetchall(sql_query,)[0][0]

                sql_id = 'select identity from orth_%s.%s where (locus_a="%s" and locus_b="%s") or (locus_a="%s" and locus_b="%s")' % (db_name,
                                                                                                                                       group,
                                                                                                                                       locus_ref,
                                                                                                                                       locus_query,
                                                                                                                                       locus_query,
                                                                                                                                       locus_ref)
                id_value = server.adaptor.execute_and_fetchall(sql_id,)[0][0]

                identity_values.append(id_value)
            one_average_id = sum(identity_values) / float(len(identity_values))
            if not all_taxons[i_taxon] in average_id:
                average_id[all_taxons[i_taxon]] = {}
            average_id[all_taxons[i_taxon]][all_taxons[y_taxon]] = [one_average_id, len(identity_values)]
            #print "one_average_id",all_taxons[i_taxon],all_taxons[y_taxon], one_average_id, len(identity_values)

    return average_id


def get_cooccurring_groups(db_name, accession1, accession2, windows_size=10, min_number_cooccurring=4, exclude=[]):
    pass
    server, db = manipulate_biosqldb.load_db(db_name)

    reference = db.lookup(accession=accession1)
    query = db.lookup(accession=accession2)

    def get_orthogroups(record):
        groups = []
        for i in record.features:
            if i.type == 'CDS':
                try:
                    groups.append(i.qualifiers['orthogroup'][0])
                except:
                    pass

        return groups

    reference_groups = get_orthogroups(reference)
    query_groups = get_orthogroups(query)


    conserved_regions = []
    for start in range(0, len(reference_groups)):
        ref_window = reference_groups[start:start+windows_size]
        for start_query in range(0, len(query_groups)):
            query_window = query_groups[start_query:start_query+windows_size]
            #print query_window, ref_window
            identical_groups = list(set(ref_window).intersection(query_window))
            if len(identical_groups) >= min_number_cooccurring:
                for i in identical_groups:
                    if i in exclude:
                        break
                else:
                    if identical_groups not in conserved_regions:
                        conserved_regions.append(identical_groups)

    return conserved_regions

def locus_list2nucleotide_fasta(biodb,locus_list):
    from chlamdb.biosqldb import manipulate_biosqldb
    from Bio import SeqRecord
    from Bio.Seq import Seq
    from Bio import SeqIO
    server, db = manipulate_biosqldb.load_db(biodb)

    filter = '"' + '","'.join(locus_list) + '"'
    sql = 'select locus_tag, accession, start, stop, strand, organism from orthology_detail' \
          ' where locus_tag in (%s)' % (biodb, filter)



    data = server.adaptor.execute_and_fetchall(sql,)
    accession_list = list(set([i[1] for i in data]))

    accession2record = {}
    for accession in accession_list:
        #print accession
        accession2record[accession] = db.lookup(accession=accession)


    record_list = []
    for one_locus in data:
        locus_tag, accession, start, stop, strand, description = one_locus
        '''
        leng = stop-start
        seq = manipulate_biosqldb.location2sequence(server, accession, biodb, start, leng)
        if strand == -1:
            seq_obj = Seq(seq)
            seq = seq_obj.reverse_complement()
        record = SeqRecord.SeqRecord(seq,
                 id=locus_tag, name=locus_tag,
                 description=description)
        '''

        for feature in accession2record[accession].features:
            if feature.type == 'CDS':
                if feature.qualifiers['locus_tag'][0] == locus_tag:
                    seq = feature.extract(accession2record[accession].seq)
        record = SeqRecord.SeqRecord(seq,
                 id=locus_tag, name=locus_tag,
                 description=description)
        record_list.append(record)
    return record_list



def orthogroup2protein_id_list(db_name):
    server, db = manipulate_biosqldb.load_db(db_name)

    sql1 = 'select orthogroup, protein_id, locus_tag from orthology_detail;' % db_name

    data1 = server.adaptor.execute_and_fetchall(sql1)

    orthogroup2protein_id_list_dico = {}
    for i in data1:
        orthogroup = i[0]
        if i[1] == '-':
            protein_id = i[2]
        else:
            protein_id = i[1]
        if not orthogroup in orthogroup2protein_id_list_dico:
            orthogroup2protein_id_list_dico[orthogroup] = [protein_id]
        else:
            orthogroup2protein_id_list_dico[orthogroup].append(protein_id)

    return orthogroup2protein_id_list_dico








def _chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]


def update_interpro_table(biodb, locus_list, locus2ortho, i):
    server, db = manipulate_biosqldb.load_db(biodb)
    y = 0
    for locus in locus_list:
        y+=1
        #print i,y, len(locus_list), locus
        sql = 'UPDATE interpro_%s SET orthogroup="%s" WHERE  locus_tag="%s";' % (biodb, locus2ortho[locus], locus)
        server.adaptor.execute(sql,)
        server.commit()

def add_orthogroup_to_interpro_table(biodb_name):

    '''
    Drop table interpro, puis on la reconstruit
    '''

    locus2ortho = locus_tag2orthogroup(biodb_name)

    import numpy
    from multiprocessing import Process
    import time

    server, db = manipulate_biosqldb.load_db(biodb_name)
    #sql = 'ALTER TABLE interpro_%s ADD orthogroup VARCHAR(100);' % biodb_name
    sql = 'select * from interpro_%s;' % biodb_name

    interpro_data = server.adaptor.execute_and_fetchall(sql,)

    try:
        sql = 'drop table interpro_%s' % biodb_name
        server.adaptor.execute(sql,)
        server.adaptor.commit()
    except:
        pass

    sql = 'CREATE TABLE interpro_%s (accession VARCHAR(100),' \
          ' locus_tag VARCHAR(200), ' \
          ' organism VARCHAR(200),  ' \
          ' taxon_id INT,' \
          ' sequence_length INT, ' \
          ' analysis VARCHAR(100) NOT NULL, ' \
          ' signature_accession VARCHAR(100), ' \
          ' signature_description VARCHAR(1000), ' \
          ' start INT, ' \
          ' stop INT, ' \
          ' score VARCHAR(20) NOT NULL, ' \
          ' interpro_accession VARCHAR(1000) NOT NULL, ' \
          ' interpro_description VARCHAR(10000),' \
          ' GO_terms varchar(10000),' \
          ' pathways varchar(10000),' \
          ' orthogroup varchar(100))' % biodb_name
    server.adaptor.execute(sql)
    i=0
    for t, one_row in enumerate(interpro_data):
        if t % 100 == 0:
            print ("%s / %s" % (t, len(interpro_data)))
        i+=1

        orthogroup = locus2ortho[one_row[1]]


        sql = 'INSERT INTO interpro_%s(accession, locus_tag, organism, taxon_id,' \
              ' sequence_length, analysis, signature_accession, signature_description, start, ' \
              ' stop, score, interpro_accession, interpro_description, GO_terms, pathways, orthogroup) ' \
              ' values ("%s", "%s", "%s", %s, %s, "%s", "%s", "%s", %s, %s, "%s", "%s", "%s", "%s", "%s", "%s");' % (biodb_name,
                                                                                             one_row[0],
                                                                                             one_row[1],
                                                                                             one_row[2],
                                                                                             one_row[3],
                                                                                             one_row[4],
                                                                                               one_row[5],
                                                                                               one_row[6],
                                                                                               one_row[7],
                                                                                               one_row[8],
                                                                                               one_row[9],
                                                                                               one_row[10],
                                                                                               one_row[11],
                                                                                               one_row[12],
                                                                                               one_row[13],
                                                                                               one_row[14],
                                                                                               orthogroup)
        #try:

        server.adaptor.execute(sql)
        server.adaptor.commit()
        #except:
        #    pass

    '''
    locus2ortho = locus_tag2orthogroup(biodb_name)
    n_cpu = 8
    n_poc_per_list = int(numpy.ceil(len(locus2ortho)/float(n_cpu)))
    query_lists = _chunks(locus2ortho.keys(), n_poc_per_list)
    print len(query_lists)
    time.sleep(5)
    procs = []
    i = 0
    for one_list in query_lists:
        i+= 1
        proc = Process(target=update_interpro_table, args=(biodb_name, one_list, locus2ortho, i))#, out_q))
        procs.append(proc)
        proc.start()

    time.sleep(5)
    for proc in procs:
        proc.join()
    '''



    '''

    import json
    import re
    server, db = manipulate_biosqldb.load_db('chlamydia_03_15')
    taxon_id2genome_description = manipulate_biosqldb.taxon_id2genome_description(server, 'chlamydia_03_15')

    for i, accession in enumerate(taxon_id2genome_description):
        #print i, accession
        description = taxon_id2genome_description[accession]
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
        description = re.sub("Protochlamydia", "P.", description)
        description = re.sub("Chlamydia", "C.", description)
        description = re.sub("Chlamydophila", "E.", description)
        description = re.sub("Estrella", "E.", description)
        description = re.sub("Rhodopirellula", "R.", description)
        description = re.sub("Methylacidiphilum", "M.", description)
        description = re.sub(" phage", "", description)
        description = re.sub("Parachlamydia", "P.", description)
        description = re.sub("Neochlamydia", "Neo.", description)
        description = re.sub("Simkania", "S.", description)
        description = re.sub("Waddlia", "W.", description)
        description = re.sub("Pirellula", "P.", description)
        description = re.sub("Rhabdochlamydiaceae sp.", "Rhabdo", description)

        taxon_id2genome_description[accession] = description
        #taxon_id2genome_description[description] = accession[0]


    with open('genome_identity_dico.json', 'r') as f:
        genome_identity = json.load(f)
    header = '\t'
    genome_identity['316'] = {}
    for i in genome_identity:
        header += '%s\t' % taxon_id2genome_description[str(i)]

    print header
    for i in genome_identity:
        row_string = '%s\t' % taxon_id2genome_description[str(i)]
        for y in genome_identity:
            #print y,i
            try:
                row_string += '%s\t' % round(genome_identity[i][y][0],2)

            except KeyError:
                try:
                    row_string += '%s\t' % round(genome_identity[y][i][0],2)
                except KeyError:
                    row_string += '-\t'
        print row_string
'''

if __name__ == '__main__':
    #print taxonomical_form('chlamydia_03_15')
    #locus_tag2best_hit_n_taxon_ids('chlamydia_03_15', 'Rhab')
    #clean_multispecies_blastnr_record('chlamydia_03_15', create_new_sql_tables=False)

    #locus_tag2n_nr_hits('chlamydia_03_15', 'Rhab')
    #locus_tag2n_blast_bacteria('chlamydia_03_15', 'Rhab')
    calculate_average_protein_identity('kpneumo_12_15')
    #
    '''
    ribosomal_proteins = get_cooccurring_groups('chlamydia_03_15', 'Rhab', 'CP001928', 10, 9)
    merged_ribosomal_proteins = []
    for i in ribosomal_proteins:
        merged_ribosomal_proteins += i

    test = ribosomal_proteins = get_cooccurring_groups('chlamydia_03_15', 'Rhab', 'CP001928', 10, 7, merged_ribosomal_proteins)
    print len(test)
    for i in test:
        print i

    wcw_syntheny = get_cooccurring_groups('chlamydia_03_15', 'Rhab', 'CP001928',3,3)
    simkania_syntheny = get_cooccurring_groups('chlamydia_03_15', 'Rhab', 'NC_015713',3,3)
    trachomatis_syntheny = get_cooccurring_groups('chlamydia_03_15', 'Rhab', 'AE001273',3,3)
    para_proto_syntheny = get_cooccurring_groups('chlamydia_03_15', 'NC_005861', 'KNic',3,3)
    tracho_muridarum = get_cooccurring_groups('chlamydia_03_15', 'NC_002620', 'AE001273',3,3)
    print len(wcw_syntheny), len(simkania_syntheny), len(trachomatis_syntheny)
    print len(para_proto_syntheny), len(tracho_muridarum)
    '''


    import json
    import re
    server, db = manipulate_biosqldb.load_db('kpneumo_12_15')
    taxon_id2genome_description = manipulate_biosqldb.taxon_id2genome_description(server, 'kpneumo_12_15')

    for i, accession in enumerate(taxon_id2genome_description):
        #print i, accession
        description = taxon_id2genome_description[accession]
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
        description = re.sub("Protochlamydia", "P.", description)
        description = re.sub("Chlamydia", "C.", description)
        description = re.sub("Chlamydophila", "E.", description)
        description = re.sub("Estrella", "E.", description)
        description = re.sub("Rhodopirellula", "R.", description)
        description = re.sub("Methylacidiphilum", "M.", description)
        description = re.sub(" phage", "", description)
        description = re.sub("Parachlamydia", "P.", description)
        description = re.sub("Neochlamydia", "Neo.", description)
        description = re.sub("Simkania", "S.", description)
        description = re.sub("Waddlia", "W.", description)
        description = re.sub("Pirellula", "P.", description)
        description = re.sub("Rhabdochlamydiaceae sp.", "Rhabdo", description)

        taxon_id2genome_description[accession] = description
        #taxon_id2genome_description[description] = accession[0]


    with open('genome_identity_dico.json', 'r') as f:
        genome_identity = json.load(f)
    header = '\t'
    genome_identity['316'] = {}
    for i in genome_identity:
        header += '%s\t' % taxon_id2genome_description[str(i)]

    #print header
    for i in genome_identity:
        row_string = '%s\t' % taxon_id2genome_description[str(i)]
        for y in genome_identity:
            #print y,i
            try:
                row_string += '%s\t' % round(genome_identity[i][y][0],2)

            except KeyError:
                try:
                    row_string += '%s\t' % round(genome_identity[y][i][0],2)
                except KeyError:
                    row_string += '-\t'
        print (row_string)
