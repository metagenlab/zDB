#!/usr/bin/env python


class CircosAccession2biplot():

    def __init__(self, server, db,
                 biodatabase_name,
                 reference_records,
                 query_records,
                 locus_highlight,
                 out_directory):

        from chlamdb.biosqldb import manipulate_biosqldb
        from chlamdb.plots import gbk2circos
        import os
        from chlamdb.biosqldb import shell_command

        reference_accessions = []
        for record in reference_records:
            reference_accessions.append(record.id.split(".")[0])
        query_accessions = []
        for record in query_records:
            query_accessions.append(record.id.split(".")[0])

        # reference_record = gbk2circos.Record(reference_record)
        # query_record = gbk2circos.Record(query_record)

        reference_taxon_id = manipulate_biosqldb.bioentry_id2taxon_id(server,
                                                                      biodatabase_name,
                                                                      reference_accessions[0])
        query_taxon_id = manipulate_biosqldb.bioentry_id2taxon_id(server, biodatabase_name, query_accessions[0])


        #draft_contigs = gbk2circos.circos_fasta_draft("/home/trestan/Dropbox/projets/phylogenomics/results/original_21_genomes_home/classification_new_genomes/fna/Cgallinaceae.fna")


        #print record_list
        circos_files_reference = gbk2circos.orthology_circos_files(server,
                                                                   reference_records,
                                                                   reference_taxon_id,
                                                                   biodatabase_name,
                                                                   out_directory,
                                                                   locus_highlight,
                                                                   taxon_list=[reference_taxon_id,
                                                                               query_taxon_id],
                                                                   query_taxon_id=query_taxon_id,
                                                                   color_missing=False)

        circos_files_query = gbk2circos.orthology_circos_files(server,
                                                               query_records,
                                                               query_taxon_id,
                                                               biodatabase_name,
                                                               out_directory,
                                                               locus_highlight,
                                                               taxon_list=[query_taxon_id,
                                                                           reference_taxon_id],
                                                               query_taxon_id=reference_taxon_id,
                                                               color_missing=False)

        chr_spacing_list_ref = []
        if len(reference_records) > 0:
            for i in range(0, len(reference_records) - 1):
                chr_spacing_list_ref.append([reference_records[i].id, reference_records[i + 1].id])
            chr_spacing_list_ref.append([reference_records[-1].id, reference_records[0].id])

        chr_spacing_list_query = []
        if len(query_records) > 0:
            for i in range(0, len(query_records) - 1):
                chr_spacing_list_query.append([query_records[i].id, query_records[i+1].id])
            chr_spacing_list_query.append([query_records[-1].id, query_records[0].id])


        circos_reference = gbk2circos.Circos_config(circos_files_reference["contigs"],chr_spacing_list_ref)
        circos_query = gbk2circos.Circos_config(circos_files_query["contigs"], chr_spacing_list_query)

        # corcos.add_plot(circos_files["GC"], type="line", r0="1.01r", r1="1.1r", color="green", fill_color="vlgreen", thickness = "2p", z = 1, rules ="")
        # corcos.add_plot(circos_files["orthogroups"], type="line", r0="1.12r", r1= "1.22r", color="black", fill_color="red", thickness = "2p", z = 1, rules ="")

        # add plus minus genes
        circos_reference.add_highlight(circos_files_reference["plus"], fill_color="grey_a1", r1="0.98r", r0="0.95r")
        circos_reference.add_highlight(circos_files_reference["minus"], fill_color="grey_a1", r1="0.95r", r0="0.92r")

        circos_query.add_highlight(circos_files_query["plus"], fill_color="grey_a1", r1="0.98r", r0="0.95r")
        circos_query.add_highlight(circos_files_query["minus"], fill_color="grey_a1", r1="0.95r", r0="0.92r")

        config_file_reference, accessions_name_reference = self.add_gene_tracks(circos_files_reference,
                                                                                circos_reference,
                                                                                [reference_accessions[0], query_accessions[0]],
                                                                                out_directory)
        config_file_query, accessions_name_query = self.add_gene_tracks(circos_files_query,
                                                                        circos_query,
                                                                        [query_accessions[0], reference_accessions[0]],
                                                                        out_directory)





        self.reference_circos = "%s.svg" % accessions_name_reference
        self.query_circos = "%s.svg" % accessions_name_query

        config_file_reference = os.path.join(out_directory, config_file_reference)
        accessions_name_reference = accessions_name_reference

        config_file_query = os.path.join(out_directory, config_file_query)
        accessions_name_query = accessions_name_query

        (stdout, stderr, return_code) = shell_command.shell_command("circos -outputfile %s -outputdir %s -conf %s" % (accessions_name_reference, out_directory, config_file_reference))
        (stdout, stderr, return_code) = shell_command.shell_command("circos -outputfile %s -outputdir %s -conf %s" % (accessions_name_query, out_directory, config_file_query))


        #cmd1 = "inkscape -g --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose --file=%s.svg" % accessions_name_reference
        #cmd2 = "inkscape -g --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose --file=%s.svg" % accessions_name_query
        #print cmd1
        #print cmd2

        #(stdout, stderr, return_code) = shell_command.shell_command(cmd1)
        #(stdout, stderr, return_code) = shell_command.shell_command(cmd2)

        #print stdout, stderr, return_code

    # add presence/absence of orthologs
    def add_gene_tracks(self, circos_files, circos, accessions, out_dir):
        import os
        r1 = 0.90
        r0 = 0.87
        for orthofile in circos_files["genomes"]:
            #print orthofile
            circos.add_highlight(orthofile, fill_color="blue", r1="%sr" % r1, r0= "%sr" % r0)
            r1 = r1-0.05
            r0 = r0-0.05

        accessions_name = ""
        for i in accessions:
            accessions_name += "_%s" % i

        config_file = "circos_config%s.txt" % accessions_name

        config_file = os.path.join(out_dir, config_file)

        t = open(config_file, "w")
        t.write(circos.get_file())

        return (config_file, "circos" + accessions_name)



class CircosAccession2multiplot():

    def __init__(self, server, db,
                 biodatabase_name,
                 reference_records,
                 queries_accession,
                 locus_highlight,
                 out_directory,
                 draft_fasta,
                 href,
                 ordered_taxons,
                 locus2label=False,
                 show_homologs=True,
                 radius=0.75,
                 locus_highlight2=[],
                 outfile_prefix=None):

        from chlamdb.biosqldb import manipulate_biosqldb
        from chlamdb.plots import gbk2circos
        import os
        from chlamdb.biosqldb import shell_command

        reference_accessions = []
        for reference_record in reference_records:
            #print "one record", reference_record
            reference_accessions.append(reference_record.id.split(".")[0])

        if not outfile_prefix:
            outfile_prefix = ''
            for accession in reference_accessions:
                outfile_prefix += accession

        out_file = "%s.svg" % outfile_prefix
        config_file = "%s.config" % outfile_prefix
        config_file_reference = os.path.join(out_directory, config_file)

        bioentry_id2taxon_id_dict = manipulate_biosqldb.bioentry_id2taxon_id_dict(server, biodatabase_name)

        reference_taxon_id = bioentry_id2taxon_id_dict[reference_accessions[0]]

        queries_taxon_id = []
        for accession in queries_accession:
            queries_taxon_id.append(bioentry_id2taxon_id_dict[accession])

        ordered_queries_taxon_id = []
        for taxon in ordered_taxons:
            if int(taxon) in queries_taxon_id:
                ordered_queries_taxon_id.append(int(taxon))

        #print 'ordered_queries_taxon_id', ordered_queries_taxon_id
        taxon2description = manipulate_biosqldb.taxon_id2genome_description(server, biodatabase_name)
        #for i in ordered_queries_taxon_id:
        #    print taxon2description[str(i)]


        #reference_record = gbk2circos.Record(reference_record)#db.lookup(accession=reference_accession))
        #locus_highlight = []#chlamydia_core_29_33
        #for i in locus_superantigens:
        #    locus_superantigens2.append(manipulate_biosqldb.locus_tag2orthogroup_id(server, i, biodatabase_name))

        #locus_highlight = klebsiella_missing # conserved_chlamydiae
        #print record_list
        circos_files_reference = gbk2circos.orthology_circos_files(server,
                                                                   reference_records,
                                                                   reference_taxon_id,
                                                                   biodatabase_name,
                                                                   out_directory,
                                                                   locus_highlight,
                                                                   taxon_list = ordered_queries_taxon_id,
                                                                   query_taxon_id=False,
                                                                   draft_data=draft_fasta,
                                                                   draft_coordinates=False,
                                                                   locus2label=locus2label,
                                                                   show_homologs=show_homologs,
                                                                   get_orthogroup_counts=True,
                                                                   locus_highlight2=locus_highlight2)


        #print "taxon_id2description_ref", taxon_id2description_reference
        chr_spacing_list = []
        #print "reference_records", len(reference_records), reference_records
        if len(reference_records) == 1 and draft_fasta[0] is None:
            for i in range(0, len(reference_records)-1):
                chr_spacing_list.append([reference_records[i].id, reference_records[i+1].id])
            chr_spacing_list.append([reference_records[-1].id, reference_records[0].id])
        elif len(reference_records) == 2 and draft_fasta[0] is not None:

            try:
                chr_spacing_list.append([draft_fasta[0][-1][0], draft_fasta[1][0][0]])
                chr_spacing_list.append([draft_fasta[0][0][0], draft_fasta[1][-1][0]])
            except:
                chr_spacing_list.append([draft_fasta[0][-1][0], reference_records[-1].name])
                chr_spacing_list.append([draft_fasta[0][0][0], reference_records[-1].name])
        else:
            chr_spacing_list.append([reference_records[-1].id, reference_records[0].id])

        circos_reference = gbk2circos.Circos_config(circos_files_reference["contigs"], chr_spacing_list, radius=radius)

        # add plus minus genes
        circos_reference.add_highlight(circos_files_reference["plus"], fill_color="grey_a1", r1="0.98r", r0="0.95r", href=href)
        circos_reference.add_highlight(circos_files_reference["minus"], fill_color="grey_a1", r1="0.95r", r0="0.92r", href=href)

        if locus2label:
            supp = '''
                    label_snuggle             = yes

                    max_snuggle_distance            = 20r
                    show_links     = yes
                    link_dims      = 10p,88p,30p,4p,4p
                    link_thickness = 2p
                    link_color     = red

                    label_size   = 24p
                    label_font   = condensed

                    padding  = 0p
                    rpadding = 0p
                    '''

            circos_reference.add_plot(circos_files_reference["labels"], type="text", r0="1r", r1="2r", color="black", rules=supp)

        last_track = self.add_gene_tracks(circos_files_reference,
                             circos_reference,
                             queries_accession,
                             config_file_reference, href=href)


        rule = """<rule>
                condition          = var(value) < 0
                fill_color         = lred
                color = red
                </rule>

                <rule>
                condition          = var(value) > 0
                fill_color         = lblue
                color = blue
                </rule>
        """

        rule2 = """<rule>
                condition          = var(value) < 0
                fill_color         = lgreen
                color = green
                </rule>

                <rule>
                condition          = var(value) > 0
                fill_color         = lblue
                color = blue
                </rule>
        """

        conditions = circos_reference.template_rules % (rule)
        circos_reference.add_plot(circos_files_reference["GC_var"], fill_color="green", r1="%sr" % (last_track -0.02), r0= "%sr" % (last_track -0.1), type="line", rules=conditions)
        conditions = circos_reference.template_rules % (rule2)
        circos_reference.add_plot(circos_files_reference["GC_skew"], fill_color="green", r1="%sr" % (last_track -0.12), r0= "%sr" % (last_track -0.2), type="line", rules=conditions)

        t = open(config_file_reference, "w")
        t.write(circos_reference.get_file())
        t.close()
        cmd = "circos -outputfile %s -outputdir %s -conf %s" % (out_file, out_directory, config_file_reference)
        #print cmd
        (stdout, stderr, return_code) = shell_command.shell_command(cmd)

    # add presence/absence of orthologs
    def add_gene_tracks(self, circos_files, circos, accessions, out_file, href):
        import os
        r1 = 0.90 # 0.90
        r0 = 0.87 # 0.89

        for orthofile in circos_files["genomes"]:
            #print orthofile
            circos.add_plot(orthofile, type="heatmap", r1="%sr" % r1, r0= "%sr" % r0, color="blues-6-seq-3, orrd-9-seq", fill_color="", thickness = "2p", z = 1, rules ="", backgrounds="",url=href)
            #circos.add_highlight(orthofile, fill_color="ortho3", r1="%sr" % r1, r0= "%sr" % r0, href=href)
            r1 = r1-0.036 # 0.012
            r0 = r0-0.036 # 0.012

        return r0

class CircosAccession2nested_plots1():

    def __init__(self, server, db,
                 biodatabase_name,
                 reference_accession,
                 queries_accession,
                 locus_highlight,
                 out_directory,
                 draft_fasta):

        from chlamdb.biosqldb import manipulate_biosqldb
        from chlamdb.plots import gbk2circos
        import os
        from chlamdb.biosqldb import shell_command



        reference_taxon_id = manipulate_biosqldb.bioentry_id2taxon_id(server, biodatabase_name, reference_accession)


        draft_contigs = gbk2circos.circos_fasta_draft(draft_fasta)

        queries_taxon_id = []
        for taxon in queries_accession:
            queries_taxon_id.append(manipulate_biosqldb.bioentry_id2taxon_id(server, biodatabase_name, taxon))

        reference_record = gbk2circos.Record(db.lookup(accession=reference_accession))

        #print record_list
        circos_files_reference = gbk2circos.orthology_circos_files(server,
                                                                                                   [reference_record],
                                                                                                    reference_taxon_id,
                                                                                                    biodatabase_name,
                                                                                                    out_directory,
                                                                                                    locus_highlight,
                                                                                                    taxon_list = queries_taxon_id,
                                                                                                    query_taxon_id=False,
                                                                                                    draft_data=draft_contigs)

        #print "taxon_id2description_ref", taxon_id2description_reference


        circos_reference = gbk2circos.Circos_config(circos_files_reference["contigs"])

        # add plus minus genes

        out_file = "%s.svg" % reference_accession
        config_file = "%s.config" % reference_accession
        config_file_reference = os.path.join(out_directory, config_file)
        out_file = os.path.join(out_directory, out_file)

        r1_minus, r0_minus = self.add_gene_tracks(circos_files_reference,
                             circos_reference,
                             queries_accession,
                             config_file_reference,
                             circos_reference)


        #(stdout, stderr, return_code) = shell_command.shell_command("circos -outputfile %s -outputdir %s -conf %s" % (out_file, out_directory, config_file_reference))

        #cmd1 = "inkscape -g --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose --file=%s" % out_file
        #print cmd1

        #(stdout, stderr, return_code) = shell_command.shell_command(cmd1)

        #print stdout, stderr, return_code

    # add presence/absence of orthologs
    def add_gene_tracks(self, circos_files, circos, accessions, out_file, circos_obj):
        import os


        r1_plus = 0.98
        r0_plus = 0.96

        r1_minus = 0.96
        r0_minus = 0.94

        for orthofile in circos_files["genomes"]:
            #print orthofile
            r1 = r0_minus - 0.02
            r0 = r0_minus - 0.03
            circos_obj.add_highlight(circos_files["plus"], fill_color="grey_a1", r1="%sr" % r1_plus, r0="%sr" % r0_plus)
            circos_obj.add_highlight(circos_files["minus"], fill_color="grey_a1", r1="%sr" % r1_minus, r0="%sr" % r0_minus)

            circos.add_highlight(orthofile, fill_color="blue", r1="%sr" % r1, r0= "%sr" % r0)
            r1_plus -= 0.09
            r0_plus -= 0.11

            r1_minus -= 0.11
            r0_minus -= 0.13

        t = open(out_file, "w")
        t.write(circos.get_file())
        return r1_minus, r0_minus

class CircosAccession2nested_plot2():

    def __init__(self, server, db,
                 biodatabase_name,
                 reference_accession,
                 query_accession,
                 locus_highlight,
                 out_directory):

        from chlamdb.biosqldb import manipulate_biosqldb
        from chlamdb.plots import gbk2circos
        import os
        from chlamdb.biosqldb import shell_command

        reference_taxon_id = manipulate_biosqldb.bioentry_id2taxon_id(server, biodatabase_name, reference_accession)
        query_taxon_id = manipulate_biosqldb.bioentry_id2taxon_id(server, biodatabase_name, query_accession)

        reference_record = gbk2circos.Record(db.lookup(accession=reference_accession))
        query_record = gbk2circos.Record(db.lookup(accession=query_accession))

        #print record_list
        circos_files_reference = gbk2circos.orthology_circos_files(server,
                                                                                                   [reference_record],
                                                                                                    reference_taxon_id,
                                                                                                    biodatabase_name,
                                                                                                    out_directory,
                                                                                                    locus_highlight,
                                                                                                    taxon_list = [reference_taxon_id,
                                                                                                                  query_taxon_id],
                                                                                                    query_taxon_id=query_taxon_id)
        circos_files_query = gbk2circos.orthology_circos_files(server, [query_record],
                                                                                                    query_taxon_id,
                                                                                                    biodatabase_name,
                                                                                                    out_directory,
                                                                                                    locus_highlight,
                                                                                                    taxon_list=[query_taxon_id,
                                                                                                               reference_taxon_id],
                                                                                                    query_taxon_id=reference_taxon_id)

        #print "taxon_id2description_ref", taxon_id2description_reference
        #print "taxon_id2description_query", taxon_id2description_query

        circos_reference = gbk2circos.Circos_config(circos_files_reference["contigs"])

        r1_plus = 0.98
        r0_plus = 0.96

        r1_minus = 0.96
        r0_minus = 0.94

        r1 = r0_minus - 0.02
        r0 = r0_minus - 0.03

        # add plus minus genes

        circos_reference.add_highlight(circos_files_reference["plus"], fill_color="grey_a1", r1="%sr" % r1_plus, r0="%sr" % r0_plus)
        circos_reference.add_highlight(circos_files_reference["minus"], fill_color="grey_a1", r1="%sr" % r1_minus, r0="%sr" % r0_minus)

        circos_reference.add_highlight(circos_files_reference["genomes"][0], fill_color="blue", r1="%sr" % r1, r0= "%sr" % r0)

        r1_plus -= 0.09
        r0_plus -= 0.09

        r1_minus -= 0.11
        r0_minus -= 0.11

        r1 = r0_minus - 0.02
        r0 = r0_minus - 0.03

        circos_reference.add_highlight(circos_files_query["plus"], fill_color="grey_a1", r1="%sr" % r1_plus, r0="%sr" % r0_plus)
        circos_reference.add_highlight(circos_files_query["minus"], fill_color="grey_a1", r1="%sr" % r1_minus, r0="%sr" % r0_minus)

        circos_reference.add_highlight(circos_files_query["genomes"][0], fill_color="blue", r1="%sr" % r1, r0= "%sr" % r0)


        out_file = "%s.svg" % reference_accession
        config_file = "%s.config" % reference_accession
        config_file_reference = os.path.join(out_directory, config_file)
        #out_file = os.path.join(out_directory, out_file)

        #print "out file", config_file_reference

        t = open(config_file_reference, "w")
        t.write(circos_reference.get_file())


        (stdout, stderr, return_code) = shell_command.shell_command("circos -outputfile %s -outputdir %s -conf %s" % (out_file, out_directory, config_file_reference))

        #cmd1 = "inkscape -g --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose --file=%s" % out_file
        #print cmd1

        #(stdout, stderr, return_code) = shell_command.shell_command(cmd1)

        #print stdout, stderr, return_code


class CircosAccession2blastnr_plot():

    def __init__(self, server,
                 biodatabase_name,
                 reference_records,
                 out_directory,
                 locus_highlight=[],
                 queries_accession=[],
                 exclude_family=False,
                 taxon_list=False,
                 highlight_BBH=False):

        from chlamdb.biosqldb import manipulate_biosqldb
        from chlamdb.plots import gbk2circos
        import os
        from chlamdb.biosqldb import shell_command

        reference_accessions = []
        #print reference_records
        for reference_record in reference_records:
            #print "one record", reference_record
            # remove version number
            reference_accessions.append(reference_record.id.split(".")[0])

        bioentry_id2taxon_id_dict = manipulate_biosqldb.bioentry_id2taxon_id_dict(server, biodatabase_name)

        # get taxon ids
        reference_taxon_id = bioentry_id2taxon_id_dict[reference_accessions[0]]

        queries_taxon_id = []
        for accession in queries_accession:
            queries_taxon_id.append(bioentry_id2taxon_id_dict[accession])

        draft_fasta = []
        for record in reference_records:
            draft_fasta.append(gbk2circos.circos_fasta_draft_misc_features(record))
        #print '###### draft ########'
        #print draft_fasta
        #print '###### draft ########'
        #print record_list

        if highlight_BBH:
            sql = 'select locus_tag from blastnr_blastnr t1 ' \
              ' inner join biosqldb.bioentry t2 on t1.query_bioentry_id=t2.bioentry_id ' \
              ' inner join biosqldb.biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
              ' inner join blastnr_blastnr_taxonomy t4 on t1.subject_taxid=t4.taxon_id ' \
              ' inner join custom_tables_locus2seqfeature_id t5 ' \
              ' on t1.seqfeature_id=t5.seqfeature_id ' \
              ' where t1.hit_number=1 and t3.name="%s" and t4.phylum!="Chlamydiae" and t1.query_taxon_id=%s;' % (biodatabase_name,
                                                                                                             biodatabase_name,
                                                                                                             biodatabase_name,
                                                                                                             reference_taxon_id)
            BBH_color = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
        else:
            BBH_color  = []
        #print 'BBH_COLOR!!', len(BBH_color)
        circos_files_reference = gbk2circos.orthology_circos_files(server,
                                                                   reference_records,
                                                                   reference_taxon_id,
                                                                   biodatabase_name,
                                                                   out_directory,
                                                                   locus_highlight=BBH_color,
                                                                   taxon_list = queries_taxon_id,
                                                                   query_taxon_id=False,
                                                                   draft_data=draft_fasta)


        #print '################ draft fasta ################'
        #print draft_fasta

        if len(draft_fasta) == 2:
            if draft_fasta[0] is None and draft_fasta[1] is None:
                draft_fasta = False

        # add spacing between chromosome and plasmids
        chr_spacing_list = []
        #print "reference_records", len(reference_records), reference_records
        if len(reference_records) > 0 and draft_fasta is False:
            for i in range(0, len(reference_records)-1):
                chr_spacing_list.append([reference_records[i].id, reference_records[i+1].id])
            chr_spacing_list.append([reference_records[-1].id, reference_records[0].id])

        #print reference_records[-1].name
        #print draft_fasta[0][-1][0]
        elif len(reference_records) == 2 and draft_fasta is not False:
            try:
                chr_spacing_list.append([draft_fasta[0][-1][0], draft_fasta[1][0][0]])
                chr_spacing_list.append([draft_fasta[0][0][0], draft_fasta[1][-1][0]])
            except:
                try:
                    chr_spacing_list.append([draft_fasta[0][-1][0], reference_records[-1].name])
                    chr_spacing_list.append([draft_fasta[0][0][0], reference_records[-1].name])
                except:
                    chr_spacing_list.append([draft_fasta[1][-1][0], reference_records[0].id])
                    chr_spacing_list.append([draft_fasta[1][0][0], reference_records[0].id])
        # get circos config object
        col_file = os.path.join(out_directory, 'colors.conf')
        circos_reference = gbk2circos.Circos_config(circos_files_reference["contigs"],
                                                    chr_spacing_list,
                                                    ideogram_spacing=0,
                                                    color_files="<<include %s>>" % col_file)


        col_data = '''


violet = 197,27,138
mgrey = 217,217,217
ortho1 = 199,233,180
ortho2 = 127,205,187
ortho3 = 127,205,187
ref = 254,153,41
pred = 255,0,55
pblue = 0, 55, 255
highlight = 107, 155, 0
group_size = 187,255,167
not_conserved = 254, 29,29
chlamydiales = 24,116,205
non_chlamydiales = 0,255,255
back = 240,240,240
blue = 1,188,255
green = 27,255,1
sta2 = 255,128,0
sta1 = 54,144,192
euk = 255,131,250

        '''
        with open(col_file, 'w') as f:
            f.write(col_data)


        # add plus minus genes


        # wrinting gc files
        #gbk2circos.print_circos_GC_file(reference_records, feature_type="CDS", out_directory=out_directory, draft_data = draft_fasta)
        # writing orthogroup size files
        # writing n blast file
        # writing n blast bacteries
        # writing n blast eukaryotes
        # writing n blast archaebacteria
        # wrinting n blast chlamydiae
        # writing n blast non chlamydiae


        blastnr_files = gbk2circos.print_blasnr_circos_files(reference_records,
                                                             biodatabase_name,
                                                             out_directory,
                                                             draft_coordinates=False,
                                                             exclude_family=exclude_family,
                                                             taxon_list=taxon_list)

        '''
        all_file_names['file_n_genomes'] = os.path.join(out_directory,"circos_n_genome_presence.txt")
        all_file_names['file_n_blastnr'] = os.path.join(out_directory,"circos_n_blastnr.txt")
        all_file_names['file_n_blast_bacteria'] = os.path.join(out_directory,"circos_n_blast_bactera.txt")
        all_file_names['file_n_blast_eukaryote'] = os.path.join(out_directory,"circos_n_blast_eukaryota.txt")
        all_file_names['file_n_blast_archae'] = os.path.join(out_directory,"circos_n_blast_archae.txt")
        all_file_names['file_n_chlamydiae'] = os.path.join(out_directory,"circos_n_chlamydiae.txt")
        all_file_names['file_n_paralogs'] = os.path.join(out_directory,"circos_n_non_chlamydiae.txt")
        file_names['gc_var_file'] = out_var
        file_names['gc_skew_file'] = out_skew
        '''

        # genes
        circos_reference.add_highlight(circos_files_reference["plus"], fill_color="grey_a1", r1="0.46r", r0="0.43r")
        circos_reference.add_highlight(circos_files_reference["minus"], fill_color="grey_a1", r1="0.43r", r0="0.40r")


        ####### gc skew

        conditions = circos_reference.template_rules % (circos_reference.template_rule('var(value) < 0', 'lred') +
                                                        circos_reference.template_rule('var(value) > 0', 'lblue'))

        circos_reference.add_plot(blastnr_files['gc_skew_file'], fill_color="green", r1="0.99r", r0= "0.85r", type="line", rules=conditions, thickness=0.2)

        ####### gc var
        '''
        conditions = circos_reference.template_rules % (circos_reference.template_rule('var(value) < 0', 'lred') +
                                                        circos_reference.template_rule('var(value) > 0', 'lblue'))

        circos_reference.add_plot(blastnr_files['gc_var_file'], fill_color="green", r1="0.99r", r0= "0.85r", type="line", rules=conditions)
        '''

        ####### orthogroup size
        conditions = circos_reference.template_rules % (circos_reference.template_rule('var(value) < 25', 'not_conserved') +
                                                        circos_reference.template_rule('var(value) > 25', 'group_size'))

        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))

        circos_reference.add_plot(blastnr_files['file_n_genomes'],
                                  thickness="0.5p",
                                  fill_color="vlgreen",
                                  color="black",
                                  r1="0.55r",
                                  r0= "0.48r",
                                  type="histogram",
                                  rules=conditions,
                                  backgrounds=backgrounds)

        ####### n blastNR
        '''
        conditions = circos_reference.template_rules % (circos_reference.template_rule('var(value) < 100', 'not_conserved') +
                                                        circos_reference.template_rule('var(value) > 99', 'non_chlamydiales'))

        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))

        circos_reference.add_plot(blastnr_files['file_n_blastnr'],
                                  thickness="0.5p",
                                  fill_color="vlgreen",
                                  color="black",
                                  r1="0.60r",
                                  r0= "0.55r",
                                  type="histogram",
                                  rules=conditions,
                                  backgrounds=backgrounds)
        '''
        ####### n blast vs Eukaryotes
        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))
        circos_reference.add_plot(blastnr_files['file_n_blast_eukaryote'],
                                  thickness="0.5p",
                                  fill_color="not_conserved",
                                  color="black",
                                  r1="0.75r",
                                  r0= "0.68r",
                                  type="histogram",
                                  min="0",
                                  max="200",
                                  backgrounds=backgrounds)

        ####### n blast vs Bacteria
        '''
        conditions = circos_reference.template_rules % (circos_reference.template_rule('var(value) < 100', 'not_conserved') +
                                                        circos_reference.template_rule('var(value) > 99', 'non_chlamydiales'))
        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))
        circos_reference.add_plot(blastnr_files['file_n_blast_bacteria'],
                                  thickness="0.5p",
                                  fill_color="vlgreen",
                                  color="black",
                                  r1="0.72r",
                                  r0= "0.67r",
                                  type="histogram",
                                  rules=conditions,
                                  backgrounds=backgrounds)
        '''

        #conditions = circos_reference.template_rules % (circos_reference.template_rule("var(value) < 100", "not_conserved"),
        #                                                circos_reference.template_rule("var(value) > 99", "non_chlamydiales"))
        ####### n blast vs Bacteria
        '''
        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))
        circos_reference.add_plot(blastnr_files['file_n_blast_chlamydiae'],
                                  thickness="0.5p",
                                  fill_color="vlgreen",
                                  color="black",
                                  r1="0.78r",
                                  r0= "0.73r",
                                  type="histogram",
                                  backgrounds=backgrounds)
        '''


        #conditions = circos_reference.template_rules % (circos_reference.template_rule("var(value) < 100", "not_conserved"),
        #                                                circos_reference.template_rule("var(value) > 99", "non_chlamydiales"))


        ####### n blast vs Non Chlamydiae

        '''
        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))
        circos_reference.add_plot(blastnr_files['file_n_blast_non_chlamydiae'],
                                  thickness="0.5p",
                                  fill_color="green, blue",
                                  color="black",
                                  r1="0.84r",
                                  r0= "0.79r",
                                  type="histogram",
                                  z=1,
                                  backgrounds=backgrounds)
        '''

        ####### stacked_bacteria/chlamydiales
        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))

        #print 'stacked------------'
        circos_reference.add_plot(blastnr_files['file_stacked_chlamydiales'],
                                  thickness="0p",
                                  fill_color="sta1, sta2",
                                  r0= "0.57r",
                                  color = False,
                                  r1 = "0.66r",
                                  type="histogram",
                                  z=1,
                                  backgrounds=backgrounds)

        ####### n blast Archae
        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))

        circos_reference.add_plot(blastnr_files['file_n_blast_archae'],
                                  thickness="0.5p",
                                  fill_color="euk",
                                  r0= "0.77r",
                                  r1 = "0.85r",
                                  type="histogram",
                                  z=1,
                                  min="0",
                                  max="200",
                                  backgrounds=backgrounds)



        out_name = ''
        for accession in reference_accessions:
            out_name += accession

        out_file = "%s.svg" % out_name
        config_file = "%s.config" % out_name
        config_file_reference = os.path.join(out_directory, config_file)
        #out_file = os.path.join(out_directory, out_file)





        t = open(config_file_reference, "w")
        t.write(circos_reference.get_file())
        t.close()
        cmd = "circos -outputfile %s -outputdir %s -conf %s" % (out_file, out_directory, config_file_reference)
        #print cmd
        (stdout, stderr, return_code) = shell_command.shell_command(cmd)


if __name__ == '__main__':

    from chlamdb.biosqldb import manipulate_biosqldb


    #refernce = db.lookup(accession="AE001273") trachomatis



    #refernce =  db.lookup(accession="NC_015713") # simkania
    #plamsid = db.lookup(accession="NC_015710") # simkania plasmid
    server, db = manipulate_biosqldb.load_db('2017_06_29b_motile_chlamydiae')
    #reference = db.lookup(accession="Rht")
    #plasmid = db.lookup(accession="RhTp")

    taxon_lst = [67,1279767,1279774,1279496,48,46,55,87925,1279815,62,1279822,66,59,52,49,64,60,804807,886707,283,314,1069693,1069694,1137444,1143376,313,1172027,1172028,1035343,307,293,1279839,1279497]

    genome = db.lookup(accession="MGLZ01000000")
    #plasmid = db.lookup(accession="RhTp")

    a = CircosAccession2blastnr_plot(server,
                     '2017_06_29b_motile_chlamydiae',
                     [genome],
                     "/home/trestan/tmp/circos_test", #"/home/trestan/work/virtualshared/documents/papiers/article_rhabdo/update_09_16/papier/figures/refseq_plast",
                     taxon_list=[],
                     highlight_BBH=True)
