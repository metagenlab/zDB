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
                 href,
                 ordered_taxons):

        import manipulate_biosqldb
        import gbk2circos
        import os
        import shell_command
        import ete2

        #if draft_fasta:
        #    draft_contigs = gbk2circos.circos_fasta_draft(draft_fasta)
        #else:
        #    draft_contigs = False

        reference_accessions = []
        for reference_record in reference_records:
            #print "one record", reference_record
            reference_accessions.append(reference_record.id.split(".")[0])
        print "reference_accession", reference_accessions


        bioentry_id2taxon_id_dict = manipulate_biosqldb.bioentry_id2taxon_id_dict(server, biodatabase_name)

        reference_taxon_id = bioentry_id2taxon_id_dict[reference_accessions[0]]

        queries_taxon_id = []
        for accession in queries_accession:
            queries_taxon_id.append(bioentry_id2taxon_id_dict[accession])





        #print ordered_taxons[0], type(ordered_taxons[0]) # str
        #print queries_taxon_id[0], type(queries_taxon_id[0]) # long

        ordered_queries_taxon_id = []
        for taxon in ordered_taxons:
            if long(taxon) in queries_taxon_id:
                ordered_queries_taxon_id.append(long(taxon))
        print ordered_taxons
        print 'ordered_queries_taxon_id', ordered_queries_taxon_id







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

        klebisella_core_strict = ["group_3137",
"group_3136",
"group_3135",
"group_3134",
"group_3132",
"group_3131",
"group_3139",
"group_3138",
"group_74",
"group_3533",
"group_3532",
"group_3126",
"group_3127",
"group_3124",
"group_3125",
"group_3133",
"group_3663",
"group_3123",
"group_3128",
"group_3129",
"group_3087",
"group_3088",
"group_96",
"group_3073",
"group_159",
"group_800",
"group_83",
"group_799"]


        klebisella_core_V1 = [
'group_3126',
'group_3127',
'group_3130',
'group_3131',
'group_3132',
'group_3134',
'group_3135',
'group_3136',
'group_3137',
'group_3138',
'group_3139',
'group_340',
'group_3532',
'group_3533',
'group_3780',
'group_3781',
'group_3782',
'group_3783',
'group_3784',
'group_3785',
'group_3786',
'group_3787',
'group_3788',
'group_3789',
'group_3832',
'group_3877',
'group_3923',
'group_3924',
'group_3925',
'group_3928',
'group_3929',
'group_3930',
'group_3931',
'group_3932',
'group_3933',
'group_3934',
'group_3935',
'group_3936',
'group_3937',
'group_3938',
'group_3939',
'group_4024',
'group_4025',
'group_4026',
'group_4027',
'group_4028',
'group_4029',
'group_4044',
'group_4045',
'group_4046',
'group_4047',
'group_4048',
'group_4049',
'group_74',
'group_3124',
'group_3125',
'group_3133',
'group_3663',
'group_3742',
'group_3745',
'group_3747',
'group_3767',
'group_4050',
'group_2979',
'group_3087',
'group_3123',
'group_3128',
'group_3129',
'roup_3770',
'group_3771',
'group_3772',
'group_3773',
'group_3774',
'group_3775',
'group_3776',
'group_3777',
'group_3778',
'group_3779',
'group_3922',
'group_3926',
'group_3927',
'group_3944',
'group_3964',
'group_4030',
'group_4037',
'group_159',
'group_3073',
'group_3088',
'group_3790',
'group_3791',
'group_3844',
'group_3986',
'group_4071',
'group_4072',
'group_4073',
'group_4074',
'group_4075',
'group_800',
'group_83',
'group_96',
'group_225',
'group_2271',
'group_228',
'group_3743',
'group_3746',
'group_799',
'group_3744',
'group_3768',
'group_3769']

        klebisella_core_all_49 = [
"group_182",
"group_2413",
"group_2416",
"group_2460",
"group_3096",
"group_3097",
"group_3099",
"group_3110",
"group_3550",
"group_3554",
"group_3555",
"group_3723",
"group_3729",
"group_3821",
"group_3838",
"group_3839",
"group_3851",
"group_3852",
"group_3853",
"group_3854",
"group_3855",
"group_3856",
"group_3878",
"group_3879",
"group_3938",
"group_3939",
"group_3993",
"group_3994",
"group_3996",
"group_3997",
"group_3998",
"group_4043",
"group_4044",
"group_4045",
"group_4046",
"group_4332",
"group_4398",
"group_4399",
"group_4456",
"group_4500",
"group_4502",
"group_4503",
"group_4504",
"group_4505",
"group_4506",
"group_4507",
"group_4508",
"group_4520",
"group_4521",
"group_4522",
"group_4523",
"group_4524",
"group_4525",
"group_4526",
"group_4527",
"group_4528",
"group_4529",
"group_4540",
"group_4541",
"group_4542",
"group_4543",
"group_4544",
"group_4545",
"group_4546",
"group_4547",
"group_4548",
"group_4563",
"group_4565",
"group_4566",
"group_4636",
"group_4700",
"group_4701",
"group_4702",
"group_4703",
"group_4704",
"group_4705",
"group_4706",
"group_4707",
"group_4708",
"group_4709",
"group_4721",
"group_59",
"group_76",
"group_109",
"group_1339",
"group_2684",
"group_2686",
"group_2704",
"group_3294",
"group_3295",
"group_3296",
"group_3297",
"group_3298",
"group_3299",
"group_3546",
"group_3579",
"group_3606",
"group_3609",
"group_3620",
"group_3621",
"group_3622",
"group_3623",
"group_3624",
"group_3625",
"group_3626",
"group_3627",
"group_3628",
"group_3645",
"group_3894",
"group_4243",
"group_4269",
"group_4358",
"group_4443",
"group_4452",
"group_4453",
"group_4454",
"group_4455",
"group_4457",
"group_4458",
"group_4459",
"group_4490",
"group_4491",
"group_4492",
"group_4493",
"group_4494",
"group_4495",
"group_4496",
"group_4498",
"group_4499",
"group_4631",
"group_4632",
"group_4633",
"group_4634",
"group_4638",
"group_4639",
"group_4640",
"group_4645",
"group_4670",
"group_4671",
"group_4672",
"group_4674",
"group_4675",
"group_4677",
"group_4678",
"group_4679",
"group_94",
"group_107",
"group_213",
"group_235",
"group_2998",
"group_3100",
"group_3101",
"group_3102",
"group_3103",
"group_3104",
"group_3105",
"group_3106",
"group_3107",
"group_3108",
"group_3109",
"group_3300",
"group_3301",
"group_3302",
"group_3303",
"group_3304",
"group_3305",
"group_3306",
"group_3307",
"group_3735",
"group_3770",
"group_3771",
"group_3772",
"group_3774",
"group_3775",
"group_3940",
"group_3941",
"group_3942",
"group_3947",
"group_3952",
"group_4039",
"group_4116",
"group_4203",
"group_4204",
"group_4205",
"group_4207",
"group_4208",
"group_4209",
"group_4221",
"group_4342",
"group_4405",
"group_4411",
"group_4412",
"group_4414",
"group_4416",
"group_4417",
"group_4430",
"group_4431",
"group_4432",
"group_4434",
"group_4470",
"group_4471",
"group_4472",
"group_4473",
"group_4474",
"group_4475",
"group_4476",
"group_4477",
"group_4478",
"group_4530",
"group_4574",
"group_4575",
"group_4576",
"group_4577",
"group_4578",
"group_4593",
"group_4594",
"group_4595",
"group_4596",
"group_4597",
"group_4650",
"group_4651",
"group_4652",
"group_4653",
"group_4654",
"group_4655",
"group_4656",
"group_4657",
"group_4658",
"group_4659",
"group_4695",
"group_4696",
"group_4697",
"group_4699",
"group_4710",
"group_4711",
"group_4712",
"group_4713",
"group_4715",
"group_4716",
"group_4717",
"group_4718",
"group_131",
"group_1344",
"group_1345",
"group_135",
"group_172",
"group_173",
"group_192",
"group_2518",
"group_2685",
"group_2687",
"group_319",
"group_3321",
"group_3322",
"group_3323",
"group_3773",
"group_3790",
"group_3915",
"group_4059",
"group_4136",
"group_4137",
"group_4206",
"group_4283",
"group_4304",
"group_4305",
"group_4306",
"group_4307",
"group_4308",
"group_4324",
"group_4325",
"group_4361",
"group_4366",
"group_4426",
"group_4510",
"group_4511",
"group_4512",
"group_4513",
"group_4514",
"group_4515",
"group_4516",
"group_4517",
"group_4518",
"group_4519",
"group_4531",
"group_4532",
"group_4533",
"group_4534",
"group_4535",
"group_4536",
"group_4537",
"group_4538",
"group_4635",
"group_4637",
"group_4662",
"group_4663",
"group_4665",
"group_61",
"group_99",
"group_223",
"group_2714",
"group_3098",
"group_3289",
"group_3615",
"group_3616",
"group_3617",
"group_3618",
"group_3619",
"group_3632",
"group_4199",
"group_4210",
"group_4216",
"group_4270",
"group_4274",
"group_4341",
"group_4343",
"group_4344",
"group_4345",
"group_4347",
"group_4368",
"group_4369",
"group_4449",
"group_4460",
"group_4466",
"group_4467",
"group_4468",
"group_4469",
"group_4480",
"group_4481",
"group_4482",
"group_4483",
"group_4485",
"group_4487",
"group_4488",
"group_4489",
"group_4590",
"group_4591",
"group_4598",
"group_4600",
"group_4601",
"group_4604",
"group_4623",
"group_4624",
"group_4641",
"group_4642",
"group_4643",
"group_4644",
"group_4646",
"group_4648",
"group_4668",
"group_4669",
"group_4698",
"group_82",
"group_83",
"group_87",
"group_209",
"group_2688",
"group_2689",
"group_3769",
"group_3911",
"group_3912",
"group_3913",
"group_3914",
"group_4135",
"group_4253",
"group_4264",
"group_4403",
"group_4404",
"group_4429",
"group_4580",
"group_4581",
"group_4584",
"group_4588",
"group_4589",
"group_4660",
"group_4661",
"group_4664",
"group_4666",
"group_4667",
"group_4680",
"group_4681",
"group_4682",
"group_4683",
"group_4684"]

        klebisella_core = ["group_3723"
,"group_3729"
,"group_4398"
,"group_3110"
,"group_2460"
,"group_4705"
,"group_4707"
,"group_4706"
,"group_4700"
,"group_4703"
,"group_4702"
,"group_4708"
,"group_3939"
,"group_3938"
,"group_3821"
,"group_2413"
,"group_2416"
,"group_4044"
,"group_4045"
,"group_4043"
,"group_59"
,"group_4548"
,"group_4543"
,"group_4542"
,"group_4541"
,"group_4540"
,"group_4547"
,"group_4546"
,"group_4545"
,"group_4544"
,"group_4721"
,"group_3099"
,"group_3096"
,"group_3097"
,"group_76"
,"group_4529"
,"group_4528"
,"group_4525"
,"group_4524"
,"group_4527"
,"group_4526"
,"group_4521"
,"group_4520"
,"group_4523"
,"group_4522"
,"group_182"
,"group_3854"
,"group_3855"
,"group_3856"
,"group_3851"
,"group_3852"
,"group_3853"
,"group_4636"
,"group_4507"
,"group_4506"
,"group_4503"
,"group_4502"
,"group_4500"
,"group_3994"
,"group_3993"
,"group_3838"
,"group_3839"
,"group_3555"
,"group_3554"
,"group_3550"
,"group_4332"
,"group_1339"
,"group_4640"
,"group_4645"
,"group_3579"
,"group_4358"
,"group_3546"
,"group_3298"
,"group_3299"
,"group_3296"
,"group_3297"
,"group_3294"
,"group_3295"
,"group_3606"
,"group_3609"
,"group_4495"
,"group_4494"
,"group_4496"
,"group_4491"
,"group_4490"
,"group_4493"
,"group_4492"
,"group_4499"
,"group_4498"
,"group_4639"
,"group_4638"
,"group_4631"
,"group_4633"
,"group_4632"
,"group_2686"
,"group_3628"
,"group_3623"
,"group_3622"
,"group_3621"
,"group_3620"
,"group_3627"
,"group_3626"
,"group_3625"
,"group_3624"
,"group_2704"
,"group_3894"
,"group_3645"
,"group_4453"
,"group_4452"
,"group_4455"
,"group_4454"
,"group_4457"
,"group_4459"
,"group_4458"
,"group_4675"
,"group_4674"
,"group_4677"
,"group_4671"
,"group_4670"
,"group_4672"
,"group_4679"
,"group_4678"
,"group_4243"
,"group_2684"
,"group_4478"
,"group_4473"
,"group_4472"
,"group_4471"
,"group_4470"
,"group_4477"
,"group_4476"
,"group_4475"
,"group_4474"
,"group_4659"
,"group_4658"
,"group_4657"
,"group_4656"
,"group_4655"
,"group_4654"
,"group_4653"
,"group_4652"
,"group_4651"
,"group_4650"
,"group_3304"
,"group_3300"
,"group_3302"
,"group_235"
,"group_3305"
,"group_3306"
,"group_3307"
,"group_3301"
,"group_3303"
,"group_4221"
,"group_3941"
,"group_4405"
,"group_4414"
,"group_4417"
,"group_4416"
,"group_4411"
,"group_4203"
,"group_3942"
,"group_3940"
,"group_3771"
,"group_3770"
,"group_3772"
,"group_3775"
,"group_4342"
,"group_4434"
,"group_4432"
,"group_4431"
,"group_4697"
,"group_4696"
,"group_2998"
,"group_4594"
,"group_4595"
,"group_4596"
,"group_4530"
,"group_3108"
,"group_3109"
,"group_3106"
,"group_3107"
,"group_3104"
,"group_3105"
,"group_3102"
,"group_3103"
,"group_3100"
,"group_3101"
,"group_4576"
,"group_4577"
,"group_4574"
,"group_4575"
,"group_4578"
,"group_4715"
,"group_4712"
,"group_4713"
,"group_3735"
,"group_4039"
,"group_3952"
,"group_4665"
,"group_4663"
,"group_4059"
,"group_61"
,"group_4538"
,"group_4536"
,"group_4537"
,"group_4534"
,"group_4535"
,"group_4532"
,"group_4533"
,"group_4531"
,"group_4283"
,"group_1345"
,"group_192"
,"group_4510"
,"group_4511"
,"group_4512"
,"group_4513"
,"group_4514"
,"group_4515"
,"group_4516"
,"group_4517"
,"group_4518"
,"group_4519"
,"group_99"
,"group_4426"
,"group_173"
,"group_172"
,"group_3915"
,"group_2687"
,"group_4136"
,"group_4137"
,"group_4361"
,"group_4635"
,"group_4637"
,"group_135"
,"group_2685"
,"group_3790"
,"group_2518"
,"group_4325"
,"group_4324"
,"group_3322"
,"group_3323"
,"group_3321"
,"group_4216"
,"group_4698"
,"group_4347"
,"group_3632"
,"group_87"
,"group_82"
,"group_83"
,"group_4487"
,"group_4485"
,"group_4482"
,"group_4483"
,"group_4480"
,"group_4481"
,"group_4488"
,"group_4489"
,"group_4624"
,"group_3289"
,"group_4591"
,"group_3616"
,"group_3617"
,"group_3615"
,"group_3618"
,"group_3619"
,"group_4604"
,"group_4600"
,"group_4601"
,"group_2714"
,"group_4598"
,"group_4590"
,"group_4449"
,"group_4199"
,"group_4668"
,"group_4669"
,"group_3098"
,"group_4274"
,"group_4466"
,"group_4467"
,"group_4460"
,"group_4468"
,"group_4469"
,"group_4641"
,"group_4642"
,"group_4644"
,"group_4646"
,"group_4253"
,"group_4135"
,"group_3769"
,"group_4429"
,"group_4684"
,"group_4680"
,"group_4681"
,"group_4682"
,"group_4683"
,"group_2688"
,"group_2689"
,"group_3914"
,"group_3913"
,"group_4666"
,"group_4667"
,"group_4664"
,"group_4660"
,"group_4661"
,"group_4589"
,"group_4588"
,"group_4584"
,"group_4581"
,"group_4580"]


        klebsiella_real_missing = [

"group_3098"
,"group_3322"
,"group_3790"
,"group_3301"
,"group_135"
,"group_3096"
,"group_3628"
,"group_3856"
,"group_3838"
,"group_3303"
,"group_4342"
,"group_3772"
,"group_3307"
,"group_4274"
,"group_3941"
,"group_3623"
,"group_82"
,"group_3952"
,"group_3323"
,"group_4216"
,"group_3609"
,"group_3624"
,"group_3645"
,"group_3555"
,"group_3770"
,"group_3942"
,"group_1345"
,"group_3616"
,"group_4325"
,"group_3626"
,"group_3615"
,"group_3853"
,"group_3295"
,"group_3622"
,"group_3617"
,"group_3104"
,"group_182"
,"group_3621"
,"group_3855"
,"group_3821"
,"group_3103"
,"group_3851"
,"group_3099"
,"group_3294"
,"group_4135"
,"group_3771"
,"group_3839"
,"group_4039"
,"group_4059"
,"group_2684"
,"group_4136"
,"group_2687"
,"group_87"
,"group_3297"
,"group_3894"
,"group_3993"
,"group_83"
,"group_2704"
,"group_4454"
,"group_2688"
,"group_3105"
,"group_3854"
,"group_4137"
,"group_3304"
,"group_3852"
,"group_4591"
,"group_4452"
,"group_99"
,"group_4665"
,"group_4283"
,"group_3108"
,"group_3109"
,"group_3775"
,"group_3302"
,"group_4715"
,"group_3627"
,"group_3102"
,"group_2689"
,"group_3299"
,"group_2413"
,"group_4398"
,"group_4491"
,"group_61"
,"group_3306"
,"group_4324"
,"group_2685"
,"group_3100"
,"group_3625"
,"group_3723"
,"group_3550"
,"group_3939"
,"group_192"
,"group_4203"
,"group_3940"
,"group_4661"
,"group_3097"
,"group_3769"
,"group_3729"
,"group_4681"
,"group_3914"
,"group_3101"
,"group_1339"
,"group_3620"
,"group_3110"
,"group_4577"
,"group_4243"
,"group_4712"
,"group_3296"
,"group_3546"
,"group_3915"
,"group_3554"
,"group_235"
,"group_3606"
,"group_3913"
,"group_3618"
,"group_3300"
,"group_4414"
,"group_3107"
,"group_3619"
,"group_3298"
,"group_3632"
,"group_3321"
,"group_4590"
,"group_3106"]


        locus_highlight = [] # klebsiella_real_missing #klebisella_core


        locus_superantigens2 = []
        #for i in locus_superantigens:
        #    locus_superantigens2.append(manipulate_biosqldb.locus_tag2orthogroup_id(server, i, biodatabase_name))
            


        #print record_list
        circos_files_reference = gbk2circos.orthology_circos_files(server,
                                                                   reference_records,
                                                                   reference_taxon_id,
                                                                   biodatabase_name,
                                                                   out_directory,
                                                                   locus_highlight,
                                                                   taxon_list = ordered_queries_taxon_id,
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
            circos.add_plot(orthofile, type="heatmap", r1="%sr" % r1, r0= "%sr" % r0, color="blues-6-seq-3, orrd-9-seq", fill_color="", thickness = "2p", z = 1, rules ="", backgrounds="",url=href)
            #circos.add_highlight(orthofile, fill_color="ortho3", r1="%sr" % r1, r0= "%sr" % r0, href=href)
            r1 = r1-0.012
            r0 = r0-0.012

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
    server, db = manipulate_biosqldb.load_db('chlamydia_04_16')

    #refernce = db.lookup(accession="AE001273") trachomatis



    #refernce =  db.lookup(accession="NC_015713") # simkania
    #plamsid = db.lookup(accession="NC_015710") # simkania plasmid

    reference =  db.lookup(accession="Rht")


    a = CircosAccession2blastnr_plot(server,
                     'chlamydia_04_16',
                     [reference],
                     "/home/trestan/work/projets/rhabdo/circos_nr_update")
