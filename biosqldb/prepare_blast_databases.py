#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def accession2coding_density(biodb, static_dir_path):
    import manipulate_biosqldb
    import gbk2fna
    import gbk2faa
    import gbk2ffn
    import gbk2table
    import os
    from Bio import SeqIO
    import shell_command

    server, db = manipulate_biosqldb.load_db(biodb)

    sql1 = 'select distinct accession from orthology_detail_%s' % biodb

    accession_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql1,)]

    db_static_path = os.path.join(static_dir_path, biodb)
    os.mkdir(db_static_path)

    faa_path = os.path.join(db_static_path, 'faa')
    print faa_path
    os.mkdir(faa_path)

    fna_path = os.path.join(db_static_path, 'fna')
    os.mkdir(fna_path)

    ffn_path = os.path.join(db_static_path, 'ffn')
    os.mkdir(ffn_path)

    gbk_path = os.path.join(db_static_path, 'gbk')
    os.mkdir(gbk_path)

    tab_path = os.path.join(db_static_path, 'tab')
    os.mkdir(tab_path)

    for n, accession in enumerate(accession_list):
        record = db.lookup(accession=accession)
        # faa + merged

        out_name_faa = os.path.join(faa_path, accession+'.faa')
        out_name_ffn = os.path.join(ffn_path, accession+'.ffn')
        out_name_fna = os.path.join(fna_path, accession+'.fna')
        out_name_tab = os.path.join(tab_path, accession+'.tab')
        out_name_gbk = os.path.join(gbk_path, accession+'.gbk')

        gbk2faa.gbk2faa(record,
                        lformat=True,
                        outname=out_name_faa)

        # fna
        gbk2fna.gbk2fna(record, outname=out_name_fna)
        # ffn
        gbk2ffn.gbk2ffn(record,outname=out_name_ffn,locus_tag=True)
        # gbk
        with open(out_name_gbk, 'w') as f:
            SeqIO.write(record, f, 'genbank')
        # tab
        gbk2table.gbk2table(record, out_name_tab)

    # merging faa, fna and ffn
    shell_command.shell_command("cd %s; cat *faa> all.faa" % faa_path)
    shell_command.shell_command("cd %s; cat *ffn> all.ffn" % ffn_path)
    shell_command.shell_command("cd %s; cat *fna> all.fna" % fna_path)

    # formatdb
    shell_command.shell_command("cd %s; for i in `ls *faa`;do formatdb -i $i -p T; done" % faa_path)
    shell_command.shell_command("cd %s; for i in `ls *ffn`;do formatdb -i $i -p F; done" % ffn_path)
    shell_command.shell_command("cd %s; for i in `ls *fna`;do formatdb -i $i -p F; done" % fna_path)

accession2coding_density("2017_06_29b_motile_chlamydiae","/home/trestan/work/dev/django/chlamydia/assets")