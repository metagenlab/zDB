#!/usr/bin/env python


class CircosAccession2biplot():

    def __init__(self, server, db,
                 biodatabase_name,
                 reference_records,
                 query_records,
                 locus_highlight,
                 out_directory):

        import manipulate_biosqldb
        import gbk2circos
        import os
        import shell_command


        reference_accessions = []
        for record in reference_records:
            reference_accessions.append(record.id.split(".")[0])
        query_accessions = []
        for record in query_records:
            query_accessions.append(record.id.split(".")[0])

        #reference_record = gbk2circos.Record(reference_record)
        #query_record = gbk2circos.Record(query_record)

        reference_taxon_id = manipulate_biosqldb.bioentry_id2taxon_id(server, biodatabase_name, reference_accessions[0])
        query_taxon_id = manipulate_biosqldb.bioentry_id2taxon_id(server, biodatabase_name, query_accessions[0])


        #draft_contigs = gbk2circos.circos_fasta_draft("/home/trestan/Dropbox/projets/phylogenomics/results/original_21_genomes_home/classification_new_genomes/fna/Cgallinaceae.fna")


        #print record_list
        circos_files_reference = gbk2circos.orthology_circos_files(server,
                                                                                                   reference_records,
                                                                                                    reference_taxon_id,
                                                                                                    biodatabase_name,
                                                                                                    out_directory,
                                                                                                    locus_highlight,
                                                                                                    taxon_list = [reference_taxon_id,
                                                                                                                  query_taxon_id],
                                                                                                    query_taxon_id=query_taxon_id,
                                                                                                    color_missing=False)

        circos_files_query = gbk2circos.orthology_circos_files(server, query_records,
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
            for i in range(0, len(reference_records)-1):
                chr_spacing_list_ref.append([reference_records[i].id, reference_records[i+1].id])
            chr_spacing_list_ref.append([reference_records[-1].id, reference_records[0].id])

        chr_spacing_list_query = []
        if len(query_records) > 0:
            for i in range(0, len(query_records)-1):
                chr_spacing_list_query.append([query_records[i].id, query_records[i+1].id])
            chr_spacing_list_query.append([query_records[-1].id, query_records[0].id])


        circos_reference = gbk2circos.Circos_config(circos_files_reference["contigs"],chr_spacing_list_ref)
        circos_query = gbk2circos.Circos_config(circos_files_query["contigs"], chr_spacing_list_query)

        #corcos.add_plot(circos_files["GC"], type="line", r0="1.01r", r1="1.1r", color="green", fill_color="vlgreen", thickness = "2p", z = 1, rules ="")
        #corcos.add_plot(circos_files["orthogroups"], type="line", r0="1.12r", r1= "1.22r", color="black", fill_color="red", thickness = "2p", z = 1, rules ="")

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
            print orthofile
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
                 href):

        import manipulate_biosqldb
        import gbk2circos
        import os
        import shell_command

        #if draft_fasta:
        #    draft_contigs = gbk2circos.circos_fasta_draft(draft_fasta)
        #else:
        #    draft_contigs = False

        reference_accessions = []
        for reference_record in reference_records:
            print "one record", reference_record
            reference_accessions.append(reference_record.id.split(".")[0])
        print "reference_accession", reference_accessions


        bioentry_id2taxon_id_dict = manipulate_biosqldb.bioentry_id2taxon_id_dict(server, biodatabase_name)

        reference_taxon_id = bioentry_id2taxon_id_dict[reference_accessions[0]]

        queries_taxon_id = []
        for accession in queries_accession:
            queries_taxon_id.append(bioentry_id2taxon_id_dict[accession])

        #reference_record = gbk2circos.Record(reference_record)#db.lookup(accession=reference_accession))

        locus_superantigens = ["SaC_00371",
        "SaC_00397",
        "SaC_00398",
        "SaC_00399",
        "SaC_00400",
        "SaC_00401",
        "SaC_00402",
        "SaC_00403",
        "SaC_00404",
        "SaC_00405",
        "SaC_00406",
        "SaC_00408",
        "SaC_01052",
        "SaC_01053",
        "SaC_01054",
        "SaC_01555",
        "SaC_01917",
        "SaCp_00018",
        "SaCp_00020",
        "SaCp_00021"]

        locus_superantigens2 = []
        #for i in locus_superantigens:
        #    locus_superantigens2.append(manipulate_biosqldb.locus_tag2orthogroup_id(server, i, biodatabase_name))
            


        #print record_list
        circos_files_reference = gbk2circos.orthology_circos_files(server,
                                                                                                   reference_records,
                                                                                                   reference_taxon_id,
                                                                                                   biodatabase_name,
                                                                                                   out_directory,
                                                                                                   locus_superantigens2,
                                                                                                   taxon_list = queries_taxon_id,
                                                                                                   query_taxon_id=False,
                                                                                               draft_data=draft_fasta)





        #print "taxon_id2description_ref", taxon_id2description_reference
        chr_spacing_list = []
        print "reference_records", len(reference_records), reference_records
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

        print chr_spacing_list

        circos_reference = gbk2circos.Circos_config(circos_files_reference["contigs"], chr_spacing_list)

        # add plus minus genes
        circos_reference.add_highlight(circos_files_reference["plus"], fill_color="grey_a1", r1="0.98r", r0="0.95r", href=href)
        circos_reference.add_highlight(circos_files_reference["minus"], fill_color="grey_a1", r1="0.95r", r0="0.92r", href=href)

        #print "writing GC files!!!!!!!!!!!"
        #gbk2circos.print_circos_GC_file(reference_records, feature_type="CDS", out_directory=out_directory, draft_data = draft_fasta)


        out_name = ''
        for accession in reference_accessions:
            out_name += accession

        out_file = "%s.svg" % out_name
        config_file = "%s.config" % out_name
        config_file_reference = os.path.join(out_directory, config_file)
        #out_file = os.path.join(out_directory, out_file)

        self.add_gene_tracks(circos_files_reference,
                             circos_reference,
                             queries_accession,
                             config_file_reference, href=href)


        cmd = "circos -outputfile %s -outputdir %s -conf %s" % (out_file, out_directory, config_file_reference)
        print cmd
        (stdout, stderr, return_code) = shell_command.shell_command(cmd)

        #cmd1 = "inkscape -g --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose --file=%s" % out_file
        #print cmd1

        #(stdout, stderr, return_code) = shell_command.shell_command(cmd1)

        #print stdout, stderr, return_code

    # add presence/absence of orthologs
    def add_gene_tracks(self, circos_files, circos, accessions, out_file, href):
        import os
        r1 = 0.90
        r0 = 0.89
        for orthofile in circos_files["genomes"]:
            print orthofile
            circos.add_highlight(orthofile, fill_color="ortho3", r1="%sr" % r1, r0= "%sr" % r0, href=href)
            r1 = r1-0.02
            r0 = r0-0.02

        t = open(out_file, "w")
        t.write(circos.get_file())


class CircosAccession2nested_plots1():

    def __init__(self, server, db,
                 biodatabase_name,
                 reference_accession,
                 queries_accession,
                 locus_highlight,
                 out_directory,
                 draft_fasta):

        import manipulate_biosqldb
        import gbk2circos
        import os
        import shell_command



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

        self.add_gene_tracks(circos_files_reference,
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
            print orthofile
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

class CircosAccession2nested_plot2():

    def __init__(self, server, db,
                 biodatabase_name,
                 reference_accession,
                 query_accession,
                 locus_highlight,
                 out_directory):

        import manipulate_biosqldb
        import gbk2circos
        import os
        import shell_command

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

        print "out file", config_file_reference

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
                 exclude_family=False):

        import manipulate_biosqldb
        import gbk2circos
        import os
        import shell_command

        reference_accessions = []
        print reference_records
        for reference_record in reference_records:
            print "one record", reference_record
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
        print '###### draft ########'
        print draft_fasta
        print '###### draft ########'
        #print record_list
        circos_files_reference = gbk2circos.orthology_circos_files(server,
                                                                   reference_records,
                                                                   reference_taxon_id,
                                                                   biodatabase_name,
                                                                   out_directory,
                                                                   locus_highlight,
                                                                   taxon_list = queries_taxon_id,
                                                                   query_taxon_id=False,
                                                                   draft_data=draft_fasta)


        print '################ draft fasta ################'
        print draft_fasta

        # add spacing between chromosome and plasmids
        chr_spacing_list = []
        print "reference_records", len(reference_records), reference_records
        if len(reference_records) > 0 and draft_fasta is False:
            for i in range(0, len(reference_records)-1):
                chr_spacing_list.append([reference_records[i].id, reference_records[i+1].id])
            chr_spacing_list.append([reference_records[-1].id, reference_records[0].id])
        elif len(reference_records) == 2 and draft_fasta is not False:
            try:
                chr_spacing_list.append([draft_fasta[0][-1][0], draft_fasta[1][0][0]])
                chr_spacing_list.append([draft_fasta[0][0][0], draft_fasta[1][-1][0]])
            except:
                chr_spacing_list.append([draft_fasta[0][-1][0], reference_records[-1].name])
                chr_spacing_list.append([draft_fasta[0][0][0], reference_records[-1].name])

        # get circos config object
        circos_reference = gbk2circos.Circos_config(circos_files_reference["contigs"], chr_spacing_list, ideogram_spacing=4)

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


        blastnr_files = gbk2circos.print_blasnr_circos_files(reference_records, biodatabase_name, out_directory, draft_coordinates=False, exclude_family=exclude_family)

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

        circos_reference.add_highlight(circos_files_reference["plus"], fill_color="grey_a1", r1="0.48r", r0="0.45r")
        circos_reference.add_highlight(circos_files_reference["minus"], fill_color="grey_a1", r1="0.45r", r0="0.42r")



        conditions = circos_reference.template_rules % (circos_reference.template_rule('var(value) < 0', 'lred') +
                                                        circos_reference.template_rule('var(value) > 0', 'lblue'))

        circos_reference.add_plot(blastnr_files['gc_skew_file'], fill_color="green", r1="0.99r", r0= "0.9r", type="line", rules=conditions)

        conditions = circos_reference.template_rules % (circos_reference.template_rule('var(value) < 25', 'not_conserved') +
                                                        circos_reference.template_rule('var(value) > 25', 'group_size'))

        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))

        circos_reference.add_plot(blastnr_files['file_n_genomes'],
                                  thickness="0.5p",
                                  fill_color="vlgreen",
                                  color="black",
                                  r1="0.54r",
                                  r0= "0.49r",
                                  type="histogram",
                                  rules=conditions,
                                  backgrounds=backgrounds)

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


        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))
        circos_reference.add_plot(blastnr_files['file_n_blast_eukaryote'],
                                  thickness="0.5p",
                                  fill_color="not_conserved",
                                  color="black",
                                  r1="0.66r",
                                  r0= "0.61r",
                                  type="histogram",
                                  backgrounds=backgrounds)

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

        #conditions = circos_reference.template_rules % (circos_reference.template_rule("var(value) < 100", "not_conserved"),
        #                                                circos_reference.template_rule("var(value) > 99", "non_chlamydiales"))
        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))
        circos_reference.add_plot(blastnr_files['file_n_blast_chlamydiae'],
                                  thickness="0.5p",
                                  fill_color="vlgreen",
                                  color="black",
                                  r1="0.78r",
                                  r0= "0.73r",
                                  type="histogram",
                                  backgrounds=backgrounds)


        #conditions = circos_reference.template_rules % (circos_reference.template_rule("var(value) < 100", "not_conserved"),
        #                                                circos_reference.template_rule("var(value) > 99", "non_chlamydiales"))
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

        backgrounds = circos_reference.template_backgrounds % (circos_reference.template_background('back'))

        circos_reference.add_plot(blastnr_files['file_stacked_chlamydiales'],
                                  thickness="0.5p",
                                  fill_color="non_chlamydiales",
                                  color="black",
                                  r0= "0.85r",
                                  type="histogram",
                                  z=1,
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
        print cmd
        (stdout, stderr, return_code) = shell_command.shell_command(cmd)


if __name__ == '__main__':

    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db('chlamydia_03_15')

    #refernce = db.lookup(accession="AE001273") trachomatis



    #refernce =  db.lookup(accession="NC_015713") # simkania
    #plamsid = db.lookup(accession="NC_015710") # simkania plasmid

    reference =  db.lookup(accession="Rhab")


    a = CircosAccession2blastnr_plot(server,
                     'chlamydia_03_15',
                     [reference],
                     "/home/trestan/Dropbox/projets/rhabdo/circos_overview/test")