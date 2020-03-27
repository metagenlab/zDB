#!/usr/bin/python

def orthogroup2blast_hits(biodb,
                          orthogroup,
                          max_n_hits=10,
                          exclude_phylum=["Chlamydiae", 'Verrucomicrobia', 'Planctomycetes'],
                          database='blast_swissprot'):

    # get max n hits for each homolog of <orthogroup> from  other phyla
    # remove redundancy
    # return fasta file
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    exclude_filter = '"' + '","'.join(exclude_phylum) + '"'


    sql = 'select t1.locus_tag,subject_accession from orthology_detail t1 ' \
              ' inner join custom_tables_locus2seqfeature_id t2 ' \
              ' on t1.locus_tag=t2.locus_tag ' \
              ' inner join %s_%s t3 on t2.seqfeature_id=t3.seqfeature_id ' \
              ' inner join blastnr_taxonomy as t4 on t3.subject_taxid=t4.taxon_id ' \
              ' where t1.orthogroup="%s" and t4.phylum not in (%s) order by t1.locus_tag, hit_number;' % (biodb,
                                                                                                biodb,
                                                                                                database,
                                                                                                biodb,
                                                                                                orthogroup,
                                                                                                exclude_filter)

    cursor.execute(sql,)
    data = cursor.fetchall()

    locus2top_n_hits = {}

    for row in data:
        if row[0] not in locus2top_n_hits:
            locus2top_n_hits[row[0]] = [row[1]]
        else:
            if len(locus2top_n_hits[row[0]]) < max_n_hits:
                locus2top_n_hits[row[0]].append(row[1])

    return locus2top_n_hits

def swissprot_accession2fasta(accession_list):

    import urllib2
    from Bio import SeqIO
    import StringIO

    accession_list_req = ','.join(accession_list)

    link = "http://www.uniprot.org/uniprot/?query=%s&format=fasta" % (accession_list_req)

    try:
        req = urllib2.Request(link)
        page = urllib2.urlopen(req)
    except urllib2.HTTPError as e:
        print (accession_list_req)
        import time
        print ('connexion error, trying again')
        print (accession_list)
        time.sleep(60)
        return swissprot_accession2fasta(accession_list)

    handle = StringIO.StringIO(page.read().decode('utf-8'))

    records = [i for i in SeqIO.parse(handle, 'fasta')]

    return records

def refseq_accession2fasta(accession_list):
    from Bio import Entrez
    from Bio import SeqIO
    import urllib2

    Entrez.email = "trestan.pillonel@unil.ch"
    try:
        handle = Entrez.efetch(db='protein', id=','.join(accession_list), rettype="fasta", retmode="text")
    except urllib2.HTTPError as e:
        print (e)
        print (accession_list)
        print ('connexion problem, waiting 60s. and trying again...')
        # run again if connexion error
        import time
        time.sleep(60)
        return refseq_accession2fasta(accession_list)
    except urllib2.URLError:
        print (e)
        print (accession_list)
        print ('connexion problem, waiting 60s. and trying again...')
        # run again if connexion error
        import time
        time.sleep(60)
        return refseq_accession2fasta(accession_list)

    records = [i for i in SeqIO.parse(handle, "fasta")]

    return records

def orthogroup2locus_and_sequences(biodb,orthogroup):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, translation from orthology_detail where orthogroup="%s"' % (biodb, orthogroup)

    data = server.adaptor.execute_and_fetchall(sql,)

    return data


def orthogroup2alignment_closest(orthogroup,
                                 biodb,
                                 max_n_hits_uniprot=2,
                                 max_n_hits_refseq=2,
                                 exclude_phylum=["Chlamydiae", 'Verrucomicrobia', 'Planctomycetes', 'Lentisphaerae'],
                                 outname="out.fa",
                                 swissprot=True,
                                 refseq=True):
    from promer2circos import chunks

    if swissprot:
        print ('getting swissprot best hits...')
        # get uniprot record of the <max_n_hits_uniprot> best hits of each <orthogroup> locus
        locus2uniprot_accession_list = orthogroup2blast_hits(biodb,
                                                             orthogroup,
                                                             max_n_hits=max_n_hits_uniprot,
                                                             exclude_phylum=exclude_phylum,
                                                             database='blast_swissprot')
        if len(locus2uniprot_accession_list) == 0:
            uniprot_sequence_records = []
            swissprot = False
        else:

            uniprot_accession_list = []
            for locus in locus2uniprot_accession_list:
                uniprot_accession_list+=locus2uniprot_accession_list[locus]

            uniprot_sequence_records = swissprot_accession2fasta(list(set(uniprot_accession_list)))

    if refseq:
        print ('getting refseq best hits...')
        # get refseq record of the <max_n_hits_refseq> best hits of each <orthogroup> locus
        locus2refseq_accession_list = orthogroup2blast_hits(biodb,
                                                            orthogroup,
                                                            max_n_hits=max_n_hits_refseq,
                                                            exclude_phylum=exclude_phylum,
                                                            database='blastnr')
        if len(locus2refseq_accession_list) == 0:
            refseq_sequence_records = []
            refseq = False
        else:
            refseq_accession_list = []
            for locus in locus2refseq_accession_list:
                refseq_accession_list+=locus2refseq_accession_list[locus]

            # download not more than 60 records at a time
            refseq_sequence_records = []
            split_lists = chunks(list(set(refseq_accession_list)), 50)
            for one_list in split_lists:
                refseq_sequence_records += refseq_accession2fasta(one_list)

    print ('writing alignment file: %s' % outname)
    ortho_sequences = orthogroup2locus_and_sequences(biodb, orthogroup)

    if len(refseq_sequence_records) > 0 or len(uniprot_sequence_records) > 0:
        with open(outname, 'w') as f:
            for locus in ortho_sequences:
                f.write(">%s\n%s\n" % (locus[0], locus[1]))
            if swissprot and len(uniprot_sequence_records) > 0:
                for record in uniprot_sequence_records:
                    name = record.name.split('|')[1]
                    f.write(">%s\n%s\n" % (name, str(record.seq)))
            if refseq and len(refseq_accession_list) > 0:
                for record in refseq_sequence_records:
                    name = record.name
                    f.write(">%s\n%s\n" % (name, str(record.seq)))
        return True
    else:
        print ('no hits for %s' % orthogroup)
        return False


def get_spaced_colors(n):

    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return ['#%02x%02x%02x' % (int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]

def plot_tree(ete3_tree,
              orthogroup,
              biodb):

    from ete3 import Tree, TreeStyle, faces, AttrFace
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    locus_list = [lf.name for lf in ete3_tree.iter_leaves()]

    filter = '"' + '","'.join(locus_list)+'"'

    print ('get uniprot taxnomy')
    sql1 = 'select subject_accession,subject_scientific_name,t2.phylum from blast_swissprot_%s t1 ' \
          ' inner join blastnr_taxonomy as t2 on t1.subject_taxid=t2.taxon_id where subject_accession in (%s);' % (biodb,
                                                                                                                   filter)
    sql1 = 'select subject_accession,subject_scientific_name,t4.phylum from orthology_detail t1 ' \
              ' inner join custom_tables_locus2seqfeature_id t2 ' \
              ' on t1.locus_tag=t2.locus_tag ' \
              ' inner join blast_swissprot_%s t3 on t2.seqfeature_id=t3.seqfeature_id ' \
              ' inner join blastnr_taxonomy as t4 on t3.subject_taxid=t4.taxon_id ' \
              ' where t1.orthogroup="%s"' % (biodb,
                                            biodb,
                                            biodb,
                                            orthogroup)
    print ('get refseq taxonomy')
    cursor.execute(sql1,)
    accession2name_and_phylum = manipulate_biosqldb.to_dict(cursor.fetchall())
    sql2 = 'select subject_accession,subject_scientific_name,t4.phylum from orthology_detail t1 ' \
              ' inner join custom_tables_locus2seqfeature_id t2 ' \
              ' on t1.locus_tag=t2.locus_tag ' \
              ' inner join blastnr_%s t3 on t2.seqfeature_id=t3.seqfeature_id ' \
              ' inner join blastnr_taxonomy as t4 on t3.subject_taxid=t4.taxon_id ' \
              ' where t1.orthogroup="%s"' % (biodb,
                                            biodb,
                                            biodb,
                                            orthogroup)

    print (sql2)
    cursor.execute(sql2,)
    accession2name_and_phylum.update(manipulate_biosqldb.to_dict(cursor.fetchall()))

    print ('plotting tree')
    phylum_list = list(set([accession2name_and_phylum[i][1] for i in accession2name_and_phylum.keys()]))

    sql = 'select locus_tag, organism from orthology_detail' % biodb
    cursor.execute(sql,)
    locus2organism = manipulate_biosqldb.to_dict(cursor.fetchall())

    phylum2col = dict(zip(phylum_list, get_spaced_colors(len(phylum_list))))

    R = ete3_tree.get_midpoint_outgroup()
    # and set it as tree outgroup
    ete3_tree.set_outgroup(R)


    for lf in ete3_tree.iter_leaves():

        try:
            swissprot_or_refseq_data = accession2name_and_phylum[lf.name]
        except:
            try:
                swissprot_or_refseq_data = accession2name_and_phylum[lf.name.split(".")[0]]
                col = phylum2col[swissprot_or_refseq_data[1]]
                lf.name = '%s|%s-%s' % (lf.name, swissprot_or_refseq_data[0], swissprot_or_refseq_data[1])

                ff = AttrFace("name", fsize=12)
                #ff.background.color = 'red'
                ff.fgcolor = col

                lf.add_face(ff, column=0)

                #nameFace = AttrFace(lf.name, fsize=30, fgcolor=phylum2col[accession2name_and_phylum[lf.name][1]])
                #faces.add_face_to_node(nameFace, lf, 0, position="branch-right")
                #
                #nameFace.border.width = 1
            except:
                col = 'red'
                try:
                    lf.name = '%s| %s' % (lf.name, locus2organism[lf.name])
                except:
                    # tryremoving version number
                    lf.name = '%s| ??' % (lf.name)
                ff = AttrFace("name", fsize=12)
                #ff.background.color = 'red'
                ff.fgcolor = col

                lf.add_face(ff, column=0)
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = True
    return ete3_tree, ts


def aafasta2phylogeny(aa_fasta, phylo=False):
    import os
    import shell_command
    import re
    from ete3 import Tree

    align_name = aa_fasta.split('.')[0] + '_mafft.fa'

    print ('aligning with mafft...')
    cmd_mafft = 'mafft --anysymbol --amino --auto --maxiterate 1000 %s > %s' % (aa_fasta, align_name)

    out, err, code = shell_command.shell_command(cmd_mafft)

    if code != 0:
        raise(err)
    if phylo:
        print ('reconstructing phylogeny with RAxML...')
        output_prefix = aa_fasta.split('.')[0]
        output_tree_name = os.path.join('RAxML_result.%s' % output_prefix)
        output_shtree_name = os.path.join('shtest_%s' % output_prefix)
        cmd_raxml = 'raxml -m PROTGAMMALG -p 12345 -s %s -n %s -c 4 -T 8;' \
              'raxml -f J -m PROTGAMMALG -s %s -p 12345 -t %s -n %s -T 8' % (align_name,
                                                                         output_prefix,
                                                                         align_name,
                                                                         output_tree_name,
                                                                         output_shtree_name)


        out, err, code = shell_command.shell_command(cmd_raxml)
        if code != 0:
            print (out, code)
            print (err)
            import sys
            sys.exit()

        nw = re.sub(":(\d+\.\d+)\[(\d+)\]", ":\\1[&&NHX:support=\\2]", open("RAxML_fastTreeSH_Support.shtest_%s" % output_prefix).read())

        t = Tree(nw, format=0)
        t.write(outfile="RAxML_fastTreeSH_Support.shtest_%s.nwk" % output_prefix, format=0)
        return t




if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    from chlamdb.biosqldb import manipulate_biosqldb
    import os
    sqlpsw = os.environ['SQLPSW']
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--orthogroup', type=str, help="orthogroup name")
    parser.add_argument("-d",'--biodb', type=str, help="biodb")
    parser.add_argument("-a",'--all_biodb', action="store_true", help="get all biodb best hits alignments")

    args = parser.parse_args()

    exclude = ['Chlamydiae'] #chlamydia_04_16: , 'Verrucomicrobia', 'Planctomycetes', 'Lentisphaerae']

    if not args.all_biodb:
        grp = str(args.orthogroup)

        alignment = orthogroup2alignment_closest(grp,
                                     args.biodb,
                                     max_n_hits_uniprot=2,
                                     max_n_hits_refseq=2,
                                     exclude_phylum=exclude,
                                     outname="%s_swiss_homologs.faa" % grp,
                                     swissprot=True,
                                     refseq=True)

        if alignment:
            t = aafasta2phylogeny("%s_swiss_homologs.faa" % grp)

            tree, ts = plot_tree(t, grp,"chlamydia_04_16")
            out_name = "%s.svg" % grp
            tree.render(out_name, tree_style=ts)
    else:
        server, db = manipulate_biosqldb.load_db(args.biodb)
        sql = 'select orthogroup, count(*) as n from orthology_detail group by orthogroup' % args.biodb

        print ('gettig orthogroup2n_hits refseq')
        orthgroup2orthogroup_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        filter = '"' + '","'.join(exclude) + '"'
        sql2 = 'select orthogroup, count(*) from ' \
               ' (select locus_tag, count(*) as n from custom_tables_locus2seqfeature_id t1 ' \
               ' inner join blastnr_blastnr as t2 on t1.seqfeature_id=t2.seqfeature_id ' \
               ' inner join blastnr_blastnr_taxonomy t3 on t2.subject_taxid=t3.taxon_id ' \
               ' where t3.phylum not in (%s) group by t1.seqfeature_id) A ' \
               ' inner join orthology_detail B on A.locus_tag=B.locus_tag ' \
               ' group by orthogroup;' % (args.biodb,
                                          args.biodb,
                                          filter,
                                          args.biodb)

        group2n_blast_refseq = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
        print ('gettig orthogroup2n_hits swissprot')
        sql3 = 'select orthogroup, count(*) from ' \
               ' (select locus_tag, count(*) as n from custom_tables_locus2seqfeature_id t1 ' \
               ' inner join blastnr_blast_swissprot as t2 on t1.seqfeature_id=t2.seqfeature_id ' \
               ' inner join blastnr_blastnr_taxonomy t3 on t2.subject_taxid=t3.taxon_id ' \
               ' where t3.phylum not in (%s) group by t1.seqfeature_id) A ' \
               ' inner join orthology_detail B on A.locus_tag=B.locus_tag ' \
               ' group by orthogroup;' % (args.biodb,
                                          args.biodb,
                                          filter,
                                          args.biodb)

        group2n_blast_swissprot = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))

        for n, group in enumerate(orthgroup2orthogroup_size):
            print ("%s / %s, %s - %s" % (n, len(orthgroup2orthogroup_size), group, orthgroup2orthogroup_size[group]))
            try:
                n_hits_swissprot = int(group2n_blast_refseq[group])
            except KeyError:
                n_hits_swissprot = 0
            try:
                n_hits_refseq = int(group2n_blast_swissprot[group])
            except KeyError:
                n_hits_refseq = 0
            if n_hits_refseq == 0 and n_hits_swissprot == 0:
                print ('not hits for %s, continue' % group)
                continue

            size = int(orthgroup2orthogroup_size[group])
            if size <= 10:
                print ('less than 10 sequences!')
                alignment = orthogroup2alignment_closest(group,
                                             args.biodb,
                                             max_n_hits_uniprot=2,
                                             max_n_hits_refseq=5,
                                             exclude_phylum=["Chlamydiae"], # , 'Verrucomicrobia', 'Planctomycetes', 'Lentisphaerae'
                                             outname="%s_swiss_homologs.faa" % group,
                                             swissprot=True,
                                             refseq=True)

            elif size <= 20:
                alignment = orthogroup2alignment_closest(group,
                                             args.biodb,
                                             max_n_hits_uniprot=2,
                                             max_n_hits_refseq=4,
                                             exclude_phylum=["Chlamydiae"], # , 'Verrucomicrobia', 'Planctomycetes', 'Lentisphaerae'
                                             outname="%s_swiss_homologs.faa" % group,
                                             swissprot=True,
                                             refseq=True)

            elif size > 100:
                alignment = orthogroup2alignment_closest(group,
                                             args.biodb,
                                             max_n_hits_uniprot=1,
                                             max_n_hits_refseq=1,
                                             exclude_phylum=["Chlamydiae"], # , 'Verrucomicrobia', 'Planctomycetes', 'Lentisphaerae'
                                             outname="%s_swiss_homologs.faa" % group,
                                             swissprot=True,
                                             refseq=True)
            else:
                alignment = orthogroup2alignment_closest(group,
                                             args.biodb,
                                             max_n_hits_uniprot=2,
                                             max_n_hits_refseq=2,
                                             exclude_phylum=["Chlamydiae"], # , 'Verrucomicrobia', 'Planctomycetes', 'Lentisphaerae'
                                             outname="%s_swiss_homologs.faa" % group,
                                             swissprot=True,
                                             refseq=True)

            #print "alignment", alignment
            #if alignment:
            #    t = aafasta2phylogeny("%s_swiss_homologs.faa" % group)
