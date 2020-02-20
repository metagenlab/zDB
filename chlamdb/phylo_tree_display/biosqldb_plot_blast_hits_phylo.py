#!/usr/bin/python


def plot_BBH_phylo(query_fasta_record, 
                   biodb,
                   asset_path, 
                   blast_type='blastp'):

    '''

    blast types:
        blastp
        blastn chromosome
        tbalen chromosome
    '''


    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    from chlamdb.biosqldb import blast_utils
    from chlamdb.phylo_tree_display import ete_motifs
    from Bio import SeqIO
    import os 
    
    try:
        label_split = True
        ordered_accession_list = [i.id.split('|')[1] for i in query_fasta_record]
    except IndexError:
        label_split = False
        ordered_accession_list = [i.id for i in query_fasta_record]
    print ("ordered_accession_list", ordered_accession_list)

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select t2.accession, t2.taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t1.name="%s" and t2.description not like "%%%%plasmid%%%%"' % biodb

    accession2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    sql = 'select locus_tag, taxon_id from orthology_detail' % biodb
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    '''
    accession2best_hits = {}
    for accession in accession2taxon_id.keys():
        print accession
        faa_path = '/home/trestan/work/dev/django/chlamydia/assets/%s/faa/%s.faa' % (biodb, accession)
        print faa_path
        one_blast = blast_utils.Blast(query_fasta_record, faa_path)
        one_blast.run_blastp()
        # keep accession and identity
        accession2best_hits[accession] = [i[0:3] for i in one_blast.best_hit_list]
    '''
    faa_path = os.path.join(asset_path, '%s/faa/all.faa' % (biodb))
    one_blast = blast_utils.Blast(query_fasta_record, faa_path)
    one_blast.run_blastp()



    locus2taxon2identity_closest = {}
    locus2taxon2locus_closest = {}

    for hit in one_blast.complete_hit_list:
        # existing locus
        if label_split:
            hit_query = hit[0].split('|')[1]
        else:
            hit_query = hit[0]
        hit_locus = hit[1]
        hit_taxon_id = str(locus_tag2taxon_id[hit_locus])


        if hit_query in locus2taxon2identity_closest:
            # existing taxon
            if hit_taxon_id in locus2taxon2identity_closest[hit_query]:
                continue
            # best hit
            else:
                locus2taxon2identity_closest[hit_query][hit_taxon_id] = float(hit[2])
                locus2taxon2locus_closest[hit_query][hit_taxon_id] = hit_locus
        else:
            # new query
            locus2taxon2identity_closest[hit_query] = {}
            locus2taxon2locus_closest[hit_query] = {}
            # best hit
            locus2taxon2identity_closest[hit_query][hit_taxon_id] = float(hit[2])
            locus2taxon2locus_closest[hit_query][hit_taxon_id] = hit_locus

    '''
    for taxon in taxon2best_hit:
        for one_hit_list in accession2best_hits[accession]:

            query_accession, hit_accession, identity = one_hit_list
            if label_split:
                label = query_accession.split('|')[1]
            else:
                label = query_accession
            if label not in locus2taxon2identity_closest:
                locus2taxon2identity_closest[label] = {}
                locus2taxon2locus_closest[label] = {}
                locus2taxon2identity_closest[label][str(accession2taxon_id[accession])] = float(identity)
                locus2taxon2locus_closest[label][str(accession2taxon_id[accession])] = hit_accession
            else:
                locus2taxon2identity_closest[label][str(accession2taxon_id[accession])] = float(identity)
                locus2taxon2locus_closest[label][str(accession2taxon_id[accession])] = hit_accession


    taxon2locus2identity_closest = {}
    query_list = []
    for accession in accession2best_hits:
        taxon2locus2identity_closest[accession2taxon_id[accession]] = {}
        for one_hit_list in accession2best_hits[accession]:
            print accession2best_hits[accession]

            query_accession, hit_accession, identity = one_hit_list
            if query_accession not in query_list:
                query_list.append(query_accession)
            taxon2locus2identity_closest[accession2taxon_id[accession]][query_accession] = identity

    print 'query list', query_list
    '''

    query_list = locus2taxon2identity_closest.keys()
    tree, style1 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                ordered_accession_list,
                                                locus2taxon2locus_closest,
                                                identity_scale=False,
                                                show_labels=True)

    tree2, style2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                ordered_accession_list,
                                                locus2taxon2identity_closest,
                                                identity_scale=True,
                                                show_labels=True)

    tree3, style3 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                ordered_accession_list,
                                                locus2taxon2identity_closest,
                                                identity_scale=True,
                                                show_labels=False)

    return tree, style1, tree2, style2, tree3, style3, locus2taxon2locus_closest



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d",'--biodb',type=str,help="biodatabase name")
    parser.add_argument("-q",'--query',type=str,help="query fasta name")

    args = parser.parse_args()

    tree, style1, tree2, style2 = plot_BBH_phylo(args.query, args.biodb)

    path = 'identity.svg' # settings.BASE_DIR + '/assets/temp/ortho_tree.svg'
    path2 = 'locus.svg'

    tree.render(path2, dpi=800, h=600, tree_style=style1)
    tree2.render(path, dpi=800, h=600, tree_style=style2)
