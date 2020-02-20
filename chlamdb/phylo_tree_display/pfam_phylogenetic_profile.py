#!/usr/bin/env python


def get_domain_taxonomy(domain_id, rank='phylum'):
    import MySQLdb
    import os
    from ete3 import NCBITaxa, Tree, TextFace,TreeStyle, StackedBarFace
    ncbi = NCBITaxa()


    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="pfam") # name of the data base
    cursor = conn.cursor()


    sql_domain_taxonomy = 'select taxid from pfam.refseq_ref_repres_genomes_domains_pfam_31 t1 ' \
                       ' inner join pfam.refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id ' \
                       ' inner join pfam.pfam_summary_version_31 t3 on t1.pfam_id=t3.hmm_id ' \
                       ' where hmm_accession="%s" group by taxid;' % (domain_id)
    cursor.execute(sql_domain_taxonomy,)

    taxid_with_domain_list = [i[0] for i in cursor.fetchall()]

    sql = 'select phylogeny from pfam.phylogeny where rank="%s"' % (rank)
    cursor.execute(sql,)
    tree_string = cursor.fetchall()[0][0]

    tree = Tree(tree_string)

    leaves_list = [i.name for i in tree.iter_leaves()]

    taxon_id2lineage = {}
    for taxon in taxid_with_domain_list:
        taxon_id2lineage[taxon] = ncbi.get_lineage(taxon)

    leaf_taxon2n_species_with_domain = {}
    for leaf_taxon in leaves_list:
        leaf_taxon2n_species_with_domain[leaf_taxon] = 0
        for taxon in taxid_with_domain_list:
            lineage = taxon_id2lineage[taxon]
            if int(leaf_taxon) in lineage:
                leaf_taxon2n_species_with_domain[leaf_taxon]+=1

                #if taxon in taxid_with_domain_list:
                #    leaf_taxon2n_species_with_domain[leaf_taxon]+=1
    return leaf_taxon2n_species_with_domain

def get_rank_summary_statistics(rank='phylum'):
    '''

    Get phylogeny from the ncbi taxonomy database given the taxon list in the table pfam.refseq_ref_repres_genomes
    Keep rank phylogeny in the table pfam.phylogeny
    Calculate genome counts for each taxon at the specified rank. Save taxid2count in the table: pfam.<rank>_leaf2n_genomes

    :param rank:
    :return:
    '''

    import MySQLdb
    import os
    from ete3 import NCBITaxa, Tree, TextFace,TreeStyle, StackedBarFace
    ncbi = NCBITaxa()


    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="pfam") # name of the data base
    cursor = conn.cursor()

    sql = 'create table if not exists pfam.phylogeny (rank varchar(400), phylogeny TEXT)'
    cursor.execute(sql,)
    conn.commit()

    sql2 = 'CREATE table if not exists pfam.leaf2n_genomes_%s(taxon_id INT, n_genomes INT)' % rank
    cursor.execute(sql2,)
    conn.commit()

    sql_taxid_list = 'select taxid from pfam.refseq_ref_repres_genomes'
    cursor.execute(sql_taxid_list,)
    taxid_list = [i[0] for i in cursor.fetchall()]

    tree = ncbi.get_topology(taxid_list, rank_limit=rank)

    taxon_id_list = [int(i.name) for i in tree.traverse("postorder")]
    taxon_id2scientific_name = ncbi.get_taxid_translator(taxon_id_list)

    sql = 'CREATE table if not exists pfam.taxid2label_%s(taxon_id INT, scientific_name TEXT, rank TEXT)' % (rank)
    cursor.execute(sql,)

    taxon_id2rank = {}
    for taxon in taxon_id2scientific_name:
        ranks = ncbi.get_rank([taxon])

        try:
            r = ranks[max(ranks.keys())]
        except:
            r = '-'
        taxon_id2rank[taxon] = r


    for taxon in taxon_id2scientific_name:
        sql = 'insert into taxid2label_%s values(%s, "%s", "%s")' % (rank,
                                                               taxon,
                                                               taxon_id2scientific_name[taxon],
                                                               taxon_id2rank[taxon])

        cursor.execute(sql,)
    conn.commit()

    collapse = ['Opisthokonta', 'Alveolata','Amoebozoa','Stramenopiles',
                'Viridiplantae','Rhodophyta', 'Trypanosomatidae', 'Viruses',
                'unclassified Bacteria', 'Leptospiraceae', 'unclassified Gammaproteobacteria',
                'unclassified Alphaproteobacteria', 'unclassified Epsilonproteobacteria',
                'unclassified Deltaproteobacteria', 'unclassified Cyanobacteria (miscellaneous)',
                 'unclassified Firmicutes sensu stricto', 'unclassified Actinobacteria (class) (miscellaneous)',
                 'unclassified Tissierellia', 'Dehalogenimonas']
    #def collapsed_leaf(node):
    #    collapse = ['Opisthokonta', 'Alveolata','Amoebozoa','Stramenopiles','Viridiplantae','Rhodophyta', 'Trypanosomatidae', 'Viruses']
    #    name = taxon_id2scientific_name[int(node.name)]
    #    if name in collapse:
    #       return True
    #    else:
    #       return False

    # colapse major euk clades some clades





    for node in tree.traverse("postorder"):
        name =  taxon_id2scientific_name[int(node.name)]
        to_detach = []
        if name in collapse:
            to_detach.extend(node.children)
            print ('ok-------------------', node.name)
        for n in to_detach:
            n.detach()
    leaves_list = [i.name for i in tree.iter_leaves()]
    leaf_taxon2n_species= {}
    leaf_taxon2n_species_with_domain = {}
    for leaf_taxon in leaves_list:
        print ('leaf', leaf_taxon)
        leaf_taxon2n_species[leaf_taxon] = 0
        leaf_taxon2n_species_with_domain[leaf_taxon] = 0
        for taxon in taxid_list:
            lineage = ncbi.get_lineage(taxon)
            if int(leaf_taxon) in lineage:
                leaf_taxon2n_species[leaf_taxon]+=1
                #if taxon in taxid_with_domain_list:
                #    leaf_taxon2n_species_with_domain[leaf_taxon]+=1
    for leaf_taxon in leaf_taxon2n_species:
        sql = 'insert into pfam.leaf2n_genomes_%s values(%s, %s)' % (rank,
                                                                           leaf_taxon,
                                                                           leaf_taxon2n_species[leaf_taxon])
        cursor.execute(sql,)
    conn.commit()

    sql = 'insert into pfam.phylogeny values("%s","%s")' % (rank, tree.write(format=1))
    cursor.execute(sql,)
    conn.commit()

def lead_reference_genome_table_into_database(genome_refseq_file=False):

    from plastnr2sqltable import insert_taxons_into_sqldb
    import MySQLdb
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="interpro") # name of the data base
    cursor = conn.cursor()

    # ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

    sql = 'create table pfam.refseq_ref_repres_genomes (assembly_id INT AUTO_INCREMENT PRIMARY KEY,' \
          ' assembly_accession varchar(200),' \
          ' bioproject varchar(200),' \
          ' biosample varchar(200),' \
          ' wgs_master varchar(200),' \
          ' refseq_category varchar(200),' \
          ' taxid INT,' \
          ' species_taxid INT,' \
          ' organism_name varchar(200),' \
          ' infraspecific_name varchar(200),' \
          ' version_status varchar(200),' \
          ' assembly_level varchar(200),' \
          ' release_type varchar(200),' \
          ' genome_rep varchar(200),' \
          ' seq_rel_date varchar(200),' \
          ' ftp_path TEXT,' \
          ' index species_taxid(species_taxid))'

    cursor.execute(sql,)
    conn.commit()

    if not genome_refseq_file:
        import urllib2
        link = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt'
        req = urllib2.Request(link)
        genome_refseq_table = urllib2.urlopen(req)
    else:
        genome_refseq_table = open(genome_refseq_file, 'r')
    taxon_update_list = []
    count = 0
    for row in genome_refseq_table:
        if row[0] == '#':
            continue
        data = row.rstrip().split('\t')
        #print data
        assembly_accession = data[0]
        bioproject = data[1]
        biosample = data[2]
        wgs_master = data[3]
        refseq_category = data[4] # reference genome / representative genome
        taxid = data[5]
        species_taxid = data[6]
        organism_name = data[7]
        infraspecific_name = data[8]
        isolate = data[9]
        version_status = data[10] # latest?
        assembly_level = data[11] # Complete genome / Chromosome / Scaffold / Contig
        release_type = data[12]
        genome_rep = data[13]
        seq_rel_date = data[14]
        ftp_path = data[19]
        if refseq_category in ['reference genome','representative genome']:
            count+=1
            print ("%s\t%s\t%s\t%s\t%s" % (assembly_accession, refseq_category, organism_name, infraspecific_name, ftp_path))
            sql = 'insert into pfam.refseq_ref_repres_genomes(assembly_accession,bioproject,biosample,wgs_master,' \
                  ' refseq_category,taxid,species_taxid,organism_name,infraspecific_name,version_status,' \
                  ' assembly_level,release_type,genome_rep,seq_rel_date,ftp_path) ' \
                  ' values("%s","%s","%s","%s","%s",%s,%s,"%s","%s",' \
                  ' "%s","%s","%s","%s","%s","%s")' % (assembly_accession,
                                                       bioproject,
                                                       biosample,
                                                       wgs_master,
                                                       refseq_category,
                                                       taxid,
                                                       species_taxid,
                                                       organism_name,
                                                       infraspecific_name,
                                                       version_status,
                                                       assembly_level,
                                                       release_type,
                                                       genome_rep,
                                                       seq_rel_date,
                                                       ftp_path)
            print (sql)
            cursor.execute(sql,)
            conn.commit()
            # check if taxon id are in blastnr_taxonomy
            sql = 'select taxon_id from blastnr_blastnr_taxonomy where taxon_id=%s'
            cursor.execute(sql % taxid,)
            try:
                taxob_id = cursor.fetchall()[0][0]
            except IndexError:
                taxon_update_list.append(taxid)

            cursor.execute(sql % species_taxid,)
            try:
                taxob_id = cursor.fetchall()[0][0]
            except IndexError:
                taxon_update_list.append(taxid)
    print ('n update taxon:', len(taxon_update_list), 'out of:',count)
    insert_taxons_into_sqldb(taxon_update_list, 300, mysql_pwd=sqlpsw)



def plot_phylum_counts(domain_id, rank='phylum',
                       colapse_low_species_counts=4,
                       remove_unlassified=True):

    '''

    1. get phylum tree
    2. foreach species => get phylum
    3. build phylum2count dictionnary
    3. plot barchart

    # merge eukaryotes into 5 main clades
    # merge virus as a single clade


    ATTENTION: no-rank groups and no-rank species...

    '''

    import MySQLdb
    import os
    from chlamdb.biosqldb import manipulate_biosqldb
    from ete3 import NCBITaxa, Tree, TextFace,TreeStyle, StackedBarFace
    ncbi = NCBITaxa()

    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="interpro") # name of the data base
    cursor = conn.cursor()

    sql = 'select * from pfam.leaf2n_genomes_%s' % rank

    cursor.execute(sql,)
    leaf_taxon2n_species = manipulate_biosqldb.to_dict(cursor.fetchall())

    leaf_taxon2n_species_with_domain = get_domain_taxonomy(domain_id, rank)

    sql = 'select phylogeny from pfam.phylogeny where rank="%s"' % (rank)

    cursor.execute(sql,)
    
    tree = Tree(cursor.fetchall()[0][0], format=1)

    sql = 'select * from pfam.taxid2label_%s' % rank
    cursor.execute(sql,)

    taxon_id2scientific_name_and_rank = manipulate_biosqldb.to_dict(cursor.fetchall())
    taxon_id2scientific_name_and_rank = {str(k):v for k,v in taxon_id2scientific_name_and_rank.items()}

    tss = TreeStyle()
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "blue"

    keep = []
    for lf in tree.iter_leaves():
        # n genomes

        if remove_unlassified:
            label = taxon_id2scientific_name_and_rank[str(lf.name)][0]
            if 'unclassified' in label:
                continue

        n_genomes = int(leaf_taxon2n_species[lf.name])
        if n_genomes > colapse_low_species_counts:
            keep.append(lf.name)
    print ('number of leaves:', len(keep))

    tree.prune(keep)

    header_list = ['Rank', 'N genomes', 'N with %s' % domain_id, 'Percentage']
    for col, header in enumerate(header_list):

        n = TextFace('%s' % (header))
        n.margin_top = 0
        n.margin_right = 1
        n.margin_left = 20
        n.margin_bottom = 1
        n.rotation = 270
        n.hz_align = 2
        n.vt_align = 2
        n.inner_background.color = "white"
        n.opacity = 1.
        tss.aligned_header.add_face(n, col)

    for lf in tree.iter_leaves():
        # n genomes

        n_genomes = int(leaf_taxon2n_species[lf.name])
        if n_genomes <= colapse_low_species_counts:
            continue


        n = TextFace('  %s ' % str(leaf_taxon2n_species[lf.name]))
        n.margin_top = 1
        n.margin_right = 1
        n.margin_left = 0
        n.margin_bottom = 1
        n.fsize = 7
        n.inner_background.color = "white"
        n.opacity = 1.
        lf.add_face(n, 2, position="aligned")

        # n genomes with domain
        m = TextFace('  %s ' % str(leaf_taxon2n_species_with_domain[lf.name]))
        m.margin_top = 1
        m.margin_right = 1
        m.margin_left = 0
        m.margin_bottom = 1
        m.fsize = 7
        m.inner_background.color = "white"
        m.opacity = 1.
        lf.add_face(m, 3, position="aligned")

        # rank
        ranks = ncbi.get_rank([lf.name])
        try:
            r = ranks[max(ranks.keys())]
        except:
            r = '-'
        n = TextFace('  %s ' % r, fsize=14, fgcolor='red')
        n.margin_top = 1
        n.margin_right = 1
        n.margin_left = 0
        n.margin_bottom = 1
        n.fsize = 7
        n.inner_background.color = "white"
        n.opacity = 1.
        lf.add_face(n, 1, position="aligned")

        # percent with target domain
        percentage = (float(leaf_taxon2n_species_with_domain[lf.name])/float(leaf_taxon2n_species[lf.name]))*100

        m = TextFace('  %s ' % str(round(percentage,2)))
        m.fsize = 1
        m.margin_top = 1
        m.margin_right = 1
        m.margin_left = 0
        m.margin_bottom = 1
        m.fsize = 7
        m.inner_background.color = "white"
        m.opacity = 1.
        lf.add_face(m, 4, position="aligned")


        b = StackedBarFace([percentage,
                            100-percentage],
                            width=100, height=10, colors=["#7fc97f", "white"])
        b.rotation= 0
        b.inner_border.color = "grey"
        b.inner_border.width = 0
        b.margin_right = 15
        b.margin_left = 0
        lf.add_face(b, 5, position="aligned")

        n = TextFace('%s' % taxon_id2scientific_name_and_rank[str(lf.name)][0], fgcolor = "black", fsize = 9) # , fstyle = 'italic'

        lf.name = " %s (%s)" % (taxon_id2scientific_name_and_rank[str(lf.name)][0], str(lf.name))
        n.margin_right = 10
        lf.add_face(n, 0)

    tss.show_leaf_name = False

    for node in tree.traverse("postorder"):
        try:
            r = taxon_id2scientific_name_and_rank[str(node.name)][1]
        except:
            pass
        try:
            if r in ['phylum', 'superkingdom', 'class', 'subphylum'] or taxon_id2scientific_name_and_rank[str(node.name)][0] in ['FCB group']:

                hola = TextFace("%s" % (taxon_id2scientific_name_and_rank[str(node.name)][0]))
                node.add_face(hola, column=0, position = "branch-top")
        except:
            pass
    return tree, tss
    #print tree.get_ascii(attributes=["name", "rank"])

def load_pfam_result_into_database(pfam_tab_file, pfam_version=31):
    import re
    from chlamdb.biosqldb import manipulate_biosqldb
    import MySQLdb
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="pfam") # name of the data base
    cursor = conn.cursor()


    '''
    # genome_id
    # PFAM id
    # count n distinct proteins with one or multiple time the considered domain
    # domain count
    '''
    sql = 'select hmm_accession, hmm_id from pfam.pfam_summary_version_%s' % pfam_version
    cursor.execute(sql,)
    pfam_accession2pfam_id = manipulate_biosqldb.to_dict(cursor.fetchall())
    sql = 'select assembly_id from refseq_ref_repres_genomes where assembly_accession like "%s%%%%";' % pfam_tab_file.split('.')[0]
    cursor.execute(sql,)
    assembly_id = cursor.fetchall()[0][0]


    sql = 'create table if not exists pfam.refseq_ref_repres_genomes_domains_pfam_%s (assembly_id INT, ' \
          ' pfam_id INT, ' \
          ' protein_id varchar(400),' \
          ' evalue_full FLOAT,' \
          ' score_full FLOAT,' \
          ' evalue_best FLOAT,' \
          ' score_best FLOAT,' \
          ' inc INT,' \
          ' index assembly_id(assembly_id),' \
          ' index pfam_id(pfam_id))' % pfam_version
    cursor.execute(sql,)
    conn.commit()

    with open(pfam_tab_file, 'r') as f:
        i=0
        for row in f:
            row = re.sub(' +' , '\t',row)
            i+=1
            #if i==10:
            #    import sys
            #    sys.exit()
            if row[0] == '#':
                continue
            data = row.rstrip().split('\t')

            protein_id = data[0]
            pfam_hmm_accession = data[3]
            evalue_full = data[4]
            score_full = data[5]
            evalue_best = data[7]
            score_best = data[8]
            inc = data[16]
            pfam_id = pfam_accession2pfam_id[pfam_hmm_accession]

            sql = 'insert into pfam.refseq_ref_repres_genomes_domains_pfam_%s values (%s, %s, "%s", %s, ' \
                  ' %s, %s, %s, %s)' % (pfam_version,
                                        assembly_id,
                                        pfam_id,
                                        protein_id,
                                        evalue_full,
                                        score_full,
                                        evalue_best,
                                        score_best,
                                        inc)
            cursor.execute(sql,)
        conn.commit()


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--table',type=str,help="RefSeq table")
    parser.add_argument("-p", '--plot', action='store_true', help="plot")
    parser.add_argument("-pf", '--pfam_tab_files', type=str, help="PFAM result tables (tab format)", nargs='+')
    parser.add_argument("-pf_db", '--pfam_database', type=str, help="PFAM database")

    args = parser.parse_args()
    #plot_phylum_counts()

    if args.table is not None:
        lead_reference_genome_table_into_database(genome_refseq_file=args.table)

    if args.plot is True:
        #get_domain_taxonomy("PF08486.9")
        #get_rank_summary_statistics('phylum')
        #get_rank_summary_statistics('order')
        t, tss = plot_phylum_counts("PF00646.32", 'order')
        t.render("/home/trestan/%s_%s.svg" % ("PF00646.32".split('.')[0], 'order'), tree_style=tss, w=600)
        #plot_phylum_counts("PF08486.9")
        #plot_phylum_counts("PF00023.29")
        #plot_phylum_counts("PF04564.14")

    if args.pfam_tab_files is not None:
        for n, one_table in enumerate(args.pfam_tab_files):
            print ('%s / %s -- %s' % (n, len(args.pfam_tab_files), one_table))
            load_pfam_result_into_database(one_table)

    if args.pfam_database:
        import hmm_utils

        hmm_utils.load_pfam_database_info_into_mysqldb(args.pfam_database, version=31)
