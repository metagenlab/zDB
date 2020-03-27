#!/usr/bin/env python

def insert_NOG_members_into_table(biodb,
                                  NOG_id, 
                                  NOG_list, 
                                  egggnog_version=451):

    '''

    list with eggnog identifiers of the form Taxid.SeqID


    :param NOG_list:
    :return:
    '''

    from plastnr2sqltable import insert_taxons_into_sqldb
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    # create NOG table if not exist
    sql = 'create table if not exists eggnog.NOG_members_v%s (NOG_id INT,' \
          ' taxon_id INT,' \
          ' protein_accession varchar(400),' \
          ' index NOG_id(NOG_id),' \
          ' index taxon_id(taxon_id))' % egggnog_version

    cursor.execute(sql,)
    conn.commit()

    taxon_update_list = []
    count = 0
    for protein in NOG_list:
        print ('protein', protein)
        # check if taxon_id is already into the database
        # if not add it to the list to add
        protein_taxon_id, protein_accession = protein.split('.', 1) # split on the first occurence of '.'
        sql = 'select taxon_id from blastnr_blastnr_taxonomy where taxon_id=%s'
        cursor.execute(sql % protein_taxon_id,)
        try:
            protein_taxon_id = cursor.fetchall()[0][0]
            count+=1
        except IndexError:
            taxon_update_list.append(protein_taxon_id)

        # insert the protein into the NOG table

        sql = 'insert into eggnog.NOG_members_v%s values (%s, %s, "%s")' % (egggnog_version,
                                                                          NOG_id,
                                                                          protein_taxon_id,
                                                                          protein_accession)

        cursor.execute(sql,)
    conn.commit()
    print ('n update taxon:', len(taxon_update_list), 'out of:', count)
    
    insert_taxons_into_sqldb(taxon_update_list, 300)




def load_eggnog_members_table(biodb,
                              table_file, 
                              egggnog_version=451):


    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    # ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

    sql = 'create table eggnog.NOG_table_v%s (NOG_id INT AUTO_INCREMENT PRIMARY KEY,' \
          ' group_name varchar(200),' \
          ' protein_count INT,' \
          ' species_count INT,' \
          ' COG_category varchar(100))' % egggnog_version

    cursor.execute(sql,)
    conn.commit()

    with open(table_file) as eggnog_members_table:
        for row in eggnog_members_table:
            data = row.rstrip().split('\t')
            print (data)

            #0 TaxonomicLevel|
            #1 GroupName|
            #2 ProteinCount|
            #3 SpeciesCount|
            #4 COGFunctionalCategory|
            #5 ProteinIDs

            TaxonomicLevel = data[0]
            group_name = data[1]
            protein_count = data[2]
            species_count = data[3]
            COG_category = data[4]
            ProteinIDs = data[5].split(',')

            sql = 'insert into eggnog.NOG_table_v%s(group_name, protein_count, species_count, COG_category) ' \
                  ' values("%s",%s,%s,"%s")' % (egggnog_version,
                                                group_name,
                                                protein_count,
                                                species_count,
                                                COG_category)
            print (sql)
            cursor.execute(sql,)
            conn.commit()
            sql = 'select NOG_id from eggnog.NOG_table_v%s where group_name="%s"' % (egggnog_version,
                                                                                    group_name)
            cursor.execute(sql,)
            NOG_id = cursor.fetchall()[0][0]
            insert_NOG_members_into_table(biodb,
                                          NOG_id, 
                                          ProteinIDs, 
                                          egggnog_version=egggnog_version)


def get_rank_summary_statistics(rank='phylum'):
    '''

    Get phylogeny from the ncbi taxonomy database given the taxon list in the table pfam.refseq_ref_repres_genomes
    Keep rank phylogeny in the table pfam.phylogeny
    Calculate genome counts for each taxon at the specified rank. Save taxid2count in the table: pfam.<rank>_leaf2n_genomes

    :param rank:
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb

    from ete3 import NCBITaxa, Tree, TextFace,TreeStyle, StackedBarFace
    ncbi = NCBITaxa()

    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    sql = 'create table if not exists eggnog.phylogeny (rank varchar(400), phylogeny TEXT)'
    cursor.execute(sql,)
    conn.commit()

    sql2 = 'CREATE table if not exists eggnog.leaf2n_genomes_%s(taxon_id INT, n_genomes INT)' % rank
    cursor.execute(sql2,)
    conn.commit()

    sql_taxid_list = 'select distinct taxon_id from eggnog.NOG_members_v451;'
    cursor.execute(sql_taxid_list,)
    taxid_list = [i[0] for i in cursor.fetchall()]

    tree = ncbi.get_topology(taxid_list, rank_limit=rank)

    taxon_id_list = [int(i.name) for i in tree.traverse("postorder")]
    taxon_id2scientific_name = ncbi.get_taxid_translator(taxon_id_list)

    sql = 'CREATE table if not exists eggnog.taxid2label_%s(taxon_id INT, scientific_name TEXT, rank TEXT)' % (rank)
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
        sql = 'insert into eggnog.taxid2label_%s values(%s, "%s", "%s")' % (rank,
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
        sql = 'insert into eggnog.leaf2n_genomes_%s values(%s, %s)' % (rank,
                                                                           leaf_taxon,
                                                                           leaf_taxon2n_species[leaf_taxon])
        cursor.execute(sql,)
    conn.commit()

    sql = 'insert into eggnog.phylogeny values("%s","%s")' % (rank, tree.write(format=1))
    cursor.execute(sql,)
    conn.commit()


def get_NOG_taxonomy(biodb, 
                     NOG_id, 
                     rank='phylum'):
    from ete3 import NCBITaxa, Tree, TextFace,TreeStyle, StackedBarFace
    ncbi = NCBITaxa()


    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor


    sql_domain_taxonomy = 'select distinct taxon_id from eggnog.NOG_members_v451 t1 ' \
                          ' inner join eggnog.NOG_table_v451 t2 on t1.NOG_id=t2.NOG_id ' \
                          ' where t2.group_name="%s";' % (NOG_id)
    cursor.execute(sql_domain_taxonomy,)

    taxid_with_domain_list = [i[0] for i in cursor.fetchall()]

    sql = 'select phylogeny from eggnog.phylogeny where rank="%s"' % (rank)
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


def plot_phylum_counts(biodb,
                       NOG_id,
                       rank='phylum',
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

    from chlamdb.biosqldb import manipulate_biosqldb
    from ete3 import NCBITaxa, Tree, TextFace,TreeStyle, StackedBarFace
    ncbi = NCBITaxa()

    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    sql = 'select * from eggnog.leaf2n_genomes_%s' % rank

    cursor.execute(sql,)
    leaf_taxon2n_species = manipulate_biosqldb.to_dict(cursor.fetchall())

    leaf_taxon2n_species_with_domain = get_NOG_taxonomy(biodb,
                                                        NOG_id, 
                                                        rank)

    sql = 'select phylogeny from eggnog.phylogeny where rank="%s"' % (rank)

    cursor.execute(sql,)
    tree = Tree(cursor.fetchall()[0][0], format=1)

    sql = 'select * from eggnog.taxid2label_%s' % rank
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
    print ('number of leaaves:', len(keep))

    tree.prune(keep)

    header_list = ['Rank', 'N genomes', 'N with %s' % NOG_id, 'Percentage']
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
        try:
            m = TextFace('  %s ' % str(leaf_taxon2n_species_with_domain[lf.name]))
        except:
           m = TextFace('  0 ')
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
        try:
            percentage = (float(leaf_taxon2n_species_with_domain[lf.name])/float(leaf_taxon2n_species[lf.name]))*100
        except:
            percentage = 0
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


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--table',type=str,help="eggnog members table")
    parser.add_argument("-d", '--biodb',type=str,help="biodb name")

    args = parser.parse_args()

    #load_eggnog_members_table(args.biodb,
    #                          args.table)
    #get_rank_summary_statistics(rank='order')
    
    tree, style = plot_phylum_counts(biodb,
                                     "COG4789", 
                                     rank='order',
                                     colapse_low_species_counts=0)

    tree.render("test.svg", tree_style=style)