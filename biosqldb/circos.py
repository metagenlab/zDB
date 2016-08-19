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

        ordered_queries_taxon_id = []
        for taxon in ordered_taxons:
            if long(taxon) in queries_taxon_id:
                ordered_queries_taxon_id.append(long(taxon))

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


        #locus_highlight = [] # klebsiella_real_missing #klebisella_core

        chlamydia_core_29_33 = ["group_619"
        ,"group_615"
        ,"group_616"
        ,"group_613"
        ,"group_19"
        ,"group_429"
        ,"group_428"
        ,"group_427"
        ,"group_422"
        ,"group_421"
        ,"group_420"
        ,"group_328"
        ,"group_324"
        ,"group_325"
        ,"group_326"
        ,"group_327"
        ,"group_320"
        ,"group_321"
        ,"group_322"
        ,"group_323"
        ,"group_1148"
        ,"group_128"
        ,"group_129"
        ,"group_125"
        ,"group_122"
        ,"group_123"
        ,"group_121"
        ,"group_695"
        ,"group_696"
        ,"group_1155"
        ,"group_711"
        ,"group_710"
        ,"group_715"
        ,"group_717"
        ,"group_716"
        ,"group_95"
        ,"group_94"
        ,"group_90"
        ,"group_232"
        ,"group_233"
        ,"group_839"
        ,"group_834"
        ,"group_835"
        ,"group_830"
        ,"group_831"
        ,"group_474"
        ,"group_475"
        ,"group_476"
        ,"group_477"
        ,"group_471"
        ,"group_472"
        ,"group_473"
        ,"group_500"
        ,"group_503"
        ,"group_505"
        ,"group_506"
        ,"group_507"
        ,"group_508"
        ,"group_509"
        ,"group_286"
        ,"group_285"
        ,"group_283"
        ,"group_975"
        ,"group_1026"
        ,"group_1024"
        ,"group_1023"
        ,"group_1021"
        ,"group_1028"
        ,"group_231"
        ,"group_218"
        ,"group_219"
        ,"group_581"
        ,"group_586"
        ,"group_621"
        ,"group_622"
        ,"group_299"
        ,"group_290"
        ,"group_291"
        ,"group_293"
        ,"group_295"
        ,"group_296"
        ,"group_297"
        ,"group_68"
        ,"group_66"
        ,"group_64"
        ,"group_65"
        ,"group_62"
        ,"group_60"
        ,"group_61"
        ,"group_315"
        ,"group_314"
        ,"group_317"
        ,"group_316"
        ,"group_311"
        ,"group_310"
        ,"group_313"
        ,"group_312"
        ,"group_319"
        ,"group_318"
        ,"group_1113"
        ,"group_1112"
        ,"group_1110"
        ,"group_1116"
        ,"group_1115"
        ,"group_157"
        ,"group_156"
        ,"group_154"
        ,"group_153"
        ,"group_152"
        ,"group_151"
        ,"group_150"
        ,"group_479"
        ,"group_807"
        ,"group_806"
        ,"group_808"
        ,"group_399"
        ,"group_398"
        ,"group_395"
        ,"group_396"
        ,"group_391"
        ,"group_392"
        ,"group_531"
        ,"group_532"
        ,"group_535"
        ,"group_534"
        ,"group_537"
        ,"group_536"
        ,"group_277"
        ,"group_275"
        ,"group_272"
        ,"group_273"
        ,"group_271"
        ,"group_1017"
        ,"group_1015"
        ,"group_748"
        ,"group_749"
        ,"group_744"
        ,"group_745"
        ,"group_742"
        ,"group_740"
        ,"group_882"
        ,"group_228"
        ,"group_920"
        ,"group_676"
        ,"group_674"
        ,"group_675"
        ,"group_348"
        ,"group_349"
        ,"group_342"
        ,"group_343"
        ,"group_340"
        ,"group_341"
        ,"group_346"
        ,"group_347"
        ,"group_344"
        ,"group_345"
        ,"group_1122"
        ,"group_1123"
        ,"group_1126"
        ,"group_1125"
        ,"group_1128"
        ,"group_1129"
        ,"group_188"
        ,"group_189"
        ,"group_184"
        ,"group_180"
        ,"group_182"
        ,"group_183"
        ,"group_397"
        ,"group_39"
        ,"group_34"
        ,"group_36"
        ,"group_31"
        ,"group_30"
        ,"group_33"
        ,"group_32"
        ,"group_445"
        ,"group_447"
        ,"group_446"
        ,"group_441"
        ,"group_440"
        ,"group_443"
        ,"group_449"
        ,"group_859"
        ,"group_856"
        ,"group_857"
        ,"group_460"
        ,"group_104"
        ,"group_105"
        ,"group_106"
        ,"group_100"
        ,"group_101"
        ,"group_102"
        ,"group_109"
        ,"group_221"
        ,"group_223"
        ,"group_1041"
        ,"group_227"
        ,"group_735"
        ,"group_738"
        ,"group_565"
        ,"group_562"
        ,"group_563"
        ,"group_560"
        ,"group_568"
        ,"group_957"
        ,"group_279"
        ,"group_608"
        ,"group_603"
        ,"group_601"
        ,"group_607"
        ,"group_606"
        ,"group_604"
        ,"group_418"
        ,"group_419"
        ,"group_412"
        ,"group_413"
        ,"group_410"
        ,"group_411"
        ,"group_416"
        ,"group_414"
        ,"group_415"
        ,"group_96"
        ,"group_339"
        ,"group_338"
        ,"group_333"
        ,"group_332"
        ,"group_331"
        ,"group_337"
        ,"group_336"
        ,"group_335"
        ,"group_334"
        ,"group_138"
        ,"group_135"
        ,"group_134"
        ,"group_137"
        ,"group_136"
        ,"group_133"
        ,"group_132"
        ,"group_686"
        ,"group_684"
        ,"group_689"
        ,"group_688"
        ,"group_267"
        ,"group_266"
        ,"group_263"
        ,"group_262"
        ,"group_88"
        ,"group_89"
        ,"group_85"
        ,"group_81"
        ,"group_82"
        ,"group_493"
        ,"group_491"
        ,"group_495"
        ,"group_498"
        ,"group_828"
        ,"group_519"
        ,"group_518"
        ,"group_517"
        ,"group_516"
        ,"group_515"
        ,"group_514"
        ,"group_513"
        ,"group_512"
        ,"group_510"
        ,"group_1030"
        ,"group_1031"
        ,"group_1032"
        ,"group_1033"
        ,"group_1036"
        ,"group_1038"
        ,"group_1068"
        ,"group_1062"
        ,"group_763"
        ,"group_768"
        ,"group_209"
        ,"group_207"
        ,"group_206"
        ,"group_205"
        ,"group_203"
        ,"group_202"
        ,"group_201"
        ,"group_597"
        ,"group_595"
        ,"group_594"
        ,"group_593"
        ,"group_592"
        ,"group_599"
        ,"group_909"
        ,"group_904"
        ,"group_906"
        ,"group_1185"
        ,"group_652"
        ,"group_654"
        ,"group_655"
        ,"group_657"
        ,"group_970"
        ,"group_973"
        ,"group_974"
        ,"group_280"
        ,"group_59"
        ,"group_58"
        ,"group_53"
        ,"group_52"
        ,"group_57"
        ,"group_56"
        ,"group_55"
        ,"group_360"
        ,"group_361"
        ,"group_362"
        ,"group_363"
        ,"group_364"
        ,"group_365"
        ,"group_366"
        ,"group_367"
        ,"group_368"
        ,"group_369"
        ,"group_162"
        ,"group_163"
        ,"group_161"
        ,"group_164"
        ,"group_165"
        ,"group_169"
        ,"group_468"
        ,"group_463"
        ,"group_462"
        ,"group_461"
        ,"group_467"
        ,"group_466"
        ,"group_465"
        ,"group_464"
        ,"group_1029"
        ,"group_1069"
        ,"group_1061"
        ,"group_1060"
        ,"group_1065"
        ,"group_548"
        ,"group_549"
        ,"group_544"
        ,"group_545"
        ,"group_547"
        ,"group_934"
        ,"group_937"
        ,"group_932"
        ,"group_939"
        ,"group_669"
        ,"group_668"
        ,"group_1124"
        ,"group_258"
        ,"group_254"
        ,"group_255"
        ,"group_257"
        ,"group_250"
        ,"group_253"
        ,"group_494"
        ,"group_359"
        ,"group_358"
        ,"group_351"
        ,"group_350"
        ,"group_352"
        ,"group_357"
        ,"group_356"
        ,"group_193"
        ,"group_191"
        ,"group_190"
        ,"group_197"
        ,"group_196"
        ,"group_194"
        ,"group_211"
        ,"group_212"
        ,"group_22"
        ,"group_23"
        ,"group_24"
        ,"group_430"
        ,"group_431"
        ,"group_432"
        ,"group_433"
        ,"group_434"
        ,"group_435"
        ,"group_436"
        ,"group_437"
        ,"group_438"
        ,"group_439"
        ,"group_843"
        ,"group_247"
        ,"group_246"
        ,"group_1157"
        ,"group_1150"
        ,"group_1158"
        ,"group_113"
        ,"group_110"
        ,"group_117"
        ,"group_116"
        ,"group_119"
        ,"group_1074"
        ,"group_1078"
        ,"group_1079"
        ,"group_706"
        ,"group_704"
        ,"group_705"
        ,"group_709"
        ,"group_490"
        ,"group_496"
        ,"group_497"
        ,"group_448"
        ,"group_574"
        ,"group_577"
        ,"group_576"
        ,"group_573"
        ,"group_578"
        ,"group_968"
        ,"group_969"
        ,"group_962"
        ,"group_966"
        ,"group_264"
        ,"group_269"
        ,"group_268"
        ,"group_444"
        ,"group_442"
        ,"group_637"
        ,"group_632"
        ,"group_638"
        ,"group_70"
        ,"group_73"
        ,"group_74"
        ,"group_77"
        ,"group_409"
        ,"group_408"
        ,"group_401"
        ,"group_400"
        ,"group_403"
        ,"group_402"
        ,"group_405"
        ,"group_404"
        ,"group_407"
        ,"group_406"
        ,"group_394"
        ,"group_306"
        ,"group_307"
        ,"group_305"
        ,"group_303"
        ,"group_301"
        ,"group_308"
        ,"group_309"
        ,"group_393"
        ,"group_148"
        ,"group_142"
        ,"group_145"
        ,"group_147"
        ,"group_1084"
        ,"group_1086"
        ,"group_1081"
        ,"group_1080"
        ,"group_1083"
        ,"group_1088"
        ,"group_489"
        ,"group_488"
        ,"group_810"
        ,"group_558"
        ,"group_584"
        ,"group_389"
        ,"group_529"
        ,"group_520"
        ,"group_526"
        ,"group_525"
        ,"group_1000"
        ,"group_773"
        ,"group_770"
        ,"group_776"
        ,"group_893"
        ,"group_890"
        ,"group_238"
        ,"group_239"
        ,"group_481"
        ,"group_483"
        ,"group_482"
        ,"group_485"
        ,"group_484"
        ,"group_487"
        ,"group_486"
        ,"group_919"
        ,"group_913"
        ,"group_645"
        ,"group_644"
        ,"group_643"
        ,"group_642"
        ,"group_815"
        ,"group_48"
        ,"group_49"
        ,"group_42"
        ,"group_43"
        ,"group_44"
        ,"group_45"
        ,"group_993"
        ,"group_992"
        ,"group_991"
        ,"group_990"
        ,"group_997"
        ,"group_996"
        ,"group_994"
        ,"group_379"
        ,"group_378"
        ,"group_377"
        ,"group_376"
        ,"group_375"
        ,"group_374"
        ,"group_373"
        ,"group_372"
        ,"group_371"
        ,"group_370"
        ,"group_1130"
        ,"group_171"
        ,"group_170"
        ,"group_174"
        ,"group_177"
        ,"group_456"
        ,"group_457"
        ,"group_454"
        ,"group_455"
        ,"group_452"
        ,"group_453"
        ,"group_450"
        ,"group_451"
        ,"group_458"
        ,"group_459"
        ,"group_867"
        ,"group_865"
        ,"group_860"
        ,"group_384"
        ,"group_385"
        ,"group_382"
        ,"group_383"
        ,"group_380"
        ,"group_381"
        ,"group_222"
        ,"group_1076"
        ,"group_1077"
        ,"group_1071"
        ,"group_1072"
        ,"group_1073"
        ,"group_727"
        ,"group_550"
        ,"group_557"
        ,"group_788"
        ,"group_782"
        ,"group_783"
        ,"group_785"
        ,"group_2"
        ,"group_1"
        ,"group_6"
        ,"group_5"
        ,"group_4"
        ,"group_940"
        ,"group_941"
        ,"group_242"
        ,"group_241"
        ,"group_245"
        ,"group_249"
        ,"group_248"]

        #locus_highlight = []#chlamydia_core_29_33
        #for i in locus_superantigens:
        #    locus_superantigens2.append(manipulate_biosqldb.locus_tag2orthogroup_id(server, i, biodatabase_name))
            
        conserved_chlamydiae = ["RhT_00337",
"RhT_01198",
"RhT_00618",
"RhT_01318",
"RhT_00475",
"RhT_00255",
"RhT_01482",
"RhT_01328",
"RhT_00133",
"RhT_01720",
"RhT_01231",
"RhT_01623",
"RhT_01545",
"RhT_01064",
"RhT_01089",
"RhT_01700",
"RhT_00264",
"RhT_00800",
"RhT_01486",
"RhT_01218",
"RhT_00620",
"RhT_00181",
"RhT_00684",
"RhT_00072",
"RhT_00668",
"RhT_01268",
"RhT_00649",
"RhT_00101",
"RhT_00835",
"RhT_00841",
"RhT_00073",
"RhT_00660",
"RhT_01378",
"RhT_00690",
"RhT_01583",
"RhT_01054",
"RhT_00190",
"RhT_00770",
"RhT_01082",
"RhT_01484",
"RhT_00386",
"RhT_00588",
"RhT_00157",
"RhT_00675",
"RhT_01140",
"RhT_01188",
"RhT_01242",
"RhT_00929",
"RhT_00216",
"RhT_00576",
"RhT_01697",
"RhT_01271",
"RhT_00896",
"RhT_00231",
"RhT_00718",
"RhT_00334",
"RhT_00894",
"RhT_01249",
"RhT_00307",
"RhT_00766",
"RhT_01618",
"RhT_00558",
"RhT_00215",
"RhT_00148",
"RhTp_00010",
"RhT_01509",
"RhT_00049",
"RhT_01444",
"RhT_00531",
"RhT_01168",
"RhT_00128",
"RhT_01146",
"RhT_01128",
"RhT_00373",
"RhT_00776",
"RhT_00518",
"RhT_01101",
"RhT_00834",
"RhT_00917",
"RhT_00316",
"RhT_01539",
"RhT_00886",
"RhT_00403",
"RhT_01450",
"RhT_00383",
"RhT_01440",
"RhT_01243",
"RhT_00768",
"RhT_01564",
"RhT_00523",
"RhT_01232",
"RhT_00331",
"RhT_00177",
"RhT_01177",
"RhT_00527",
"RhT_00124",
"RhT_01280",
"RhT_01259",
"RhT_01353",
"RhT_00489",
"RhT_00332",
"RhT_00149",
"RhT_00671",
"RhT_00099",
"RhT_00234",
"RhT_01137",
"RhT_00214",
"RhT_00050",
"RhT_00764",
"RhT_00034",
"RhT_01363",
"RhT_00549",
"RhT_00747",
"RhT_00572",
"RhT_01600",
"RhT_00927",
"RhT_00514",
"RhT_00039",
"RhT_00333",
"RhT_00878",
"RhT_00663",
"RhT_00399",
"RhT_01451",
"RhT_00571",
"RhT_01219",
"RhT_01169",
"RhT_00393",
"RhT_00995",
"RhT_00052",
"RhT_01335",
"RhT_00293",
"RhT_01291",
"RhT_00767",
"RhT_00047",
"RhT_00996",
"RhT_00529",
"RhT_01362",
"RhT_01361",
"RhT_01009",
"RhT_00594",
"RhT_00762",
"RhT_00193",
"RhT_00294",
"RhT_00413",
"RhT_00773",
"RhT_00709",
"RhT_01441",
"RhT_00229",
"RhT_00484",
"RhT_00592",
"RhT_00411",
"RhT_01094",
"RhT_00943",
"RhT_00854",
"RhT_00581",
"RhT_01606",
"RhT_01139",
"RhT_00962",
"RhT_00813",
"RhT_00075",
"RhT_00526",
"RhT_00378",
"RhT_01125",
"RhT_00444",
"RhT_00685",
"RhT_00088",
"RhT_00253",
"RhT_01519",
"RhT_01155",
"RhT_01164",
"RhT_00691",
"RhT_00401",
"RhT_01260",
"RhT_00407",
"RhT_01602",
"RhT_00496",
"RhT_00063",
"RhT_00233",
"RhT_00228",
"RhT_00595",
"RhT_00224",
"RhT_01036",
"RhT_00178",
"RhT_00102",
"RhT_00395",
"RhT_00140",
"RhT_01124",
"RhT_01394",
"RhT_00815",
"RhT_01295",
"RhT_01236",
"RhT_01446",
"RhT_00812",
"RhT_00249",
"RhT_00400",
"RhT_00406",
"RhT_01571",
"RhT_00106",
"RhT_00584",
"RhT_01018",
"RhT_00885",
"RhT_01538",
"RhT_01652",
"RhT_01464",
"RhT_00132",
"RhT_00530",
"RhT_01447",
"RhT_01458",
"RhT_00784",
"RhT_00152",
"RhT_01059",
"RhT_00076",
"RhT_00607",
"RhT_01120",
"RhT_00416",
"RhT_01211",
"RhT_00428",
"RhT_00105",
"RhT_00636",
"RhT_01091",
"RhT_00730",
"RhT_00512",
"RhT_01033",
"RhT_00832",
"RhT_01239",
"RhT_00945",
"RhT_00631",
"RhT_01443",
"RhT_01012",
"RhT_00605",
"RhT_01480",
"RhT_00460",
"RhT_01511",
"RhT_01096",
"RhT_01184",
"RhT_01499",
"RhT_01477",
"RhT_00850",
"RhT_00573",
"RhT_00204",
"RhT_01321",
"RhT_00647",
"RhT_00577",
"RhT_00998",
"RhT_00103",
"RhT_00462",
"RhT_01313",
"RhT_01073",
"RhT_01411",
"RhT_01138",
"RhT_00305",
"RhT_01601",
"RhT_00194",
"RhT_01460",
"RhT_00402",
"RhT_00357",
"RhT_00855",
"RhT_00946",
"RhT_00982",
"RhT_00706",
"RhT_00942",
"RhT_01537",
"RhT_00404",
"RhT_00679",
"RhT_01042",
"RhT_00320",
"RhT_00084",
"RhT_00156",
"RhT_00488",
"RhT_01676",
"RhT_00129",
"RhT_01221",
"RhT_00700",
"RhT_01093",
"RhT_01266",
"RhT_00808",
"RhT_00336",
"RhT_01013",
"RhT_00814",
"RhT_00005",
"RhT_00485",
"RhT_00398",
"RhT_01233",
"RhT_01286",
"RhT_00953",
"RhT_01047",
"RhT_00318",
"RhT_01314",
"RhT_01576",
"RhT_00554",
"RhT_00557",
"RhT_01553",
"RhT_00392",
"RhT_00553",
"RhT_01290",
"RhT_00535",
"RhT_00362",
"RhT_01106",
"RhT_00254",
"RhT_01673",
"RhT_01086",
"RhT_01398",
"RhT_01148",
"RhT_00919",
"RhT_00699",
"RhT_01501",
"RhT_00874",
"RhT_01523",
"RhT_00719",
"RhT_01343",
"RhT_01105",
"RhT_01112",
"RhT_00570",
"RhT_00602",
"RhT_00191",
"RhT_00954",
"RhT_01503",
"RhT_01043",
"RhT_01234",
"RhT_01472",
"RhT_00180",
"RhT_00912",
"RhT_01463",
"RhT_00468",
"RhT_00227",
"RhT_00441",
"RhT_00350",
"RhT_01624",
"RhT_01424",
"RhT_00582",
"RhT_01302",
"RhT_00565",
"RhT_00158",
"RhT_01445",
"RhT_00997",
"RhT_00309",
"RhT_01438",
"RhT_01098",
"RhT_00068",
"RhT_01332",
"RhT_01189",
"RhT_00387",
"RhT_01665",
"RhT_00238",
"RhT_00977",
"RhT_00036",
"RhT_00723",
"RhT_00053",
"RhT_01635",
"RhT_01341",
"RhT_01418",
"RhT_01590",
"RhT_00314",
"RhT_01549",
"RhT_00733",
"RhT_00513",
"RhT_00564",
"RhT_01342",
"RhT_01037",
"RhT_00952",
"RhT_00252",
"RhT_00568",
"RhT_01427",
"RhT_01113",
"RhT_01299",
"RhT_00349",
"RhT_01584",
"RhT_01287",
"RhT_00851",
"RhT_01548",
"RhT_01052",
"RhT_01434",
"RhT_00260",
"RhT_00763",
"RhT_01465",
"RhT_00382",
"RhT_01049",
"RhT_00322",
"RhT_00409",
"RhT_01513",
"RhT_00976",
"RhT_01575",
"RhT_01309",
"RhT_00707",
"RhT_01081",
"RhT_01636",
"RhT_00285",
"RhT_01591",
"RhT_00883",
"RhT_00006",
"RhT_01634",
"RhT_01061",
"RhT_00185",
"RhT_01289",
"RhT_00887",
"RhT_00109",
"RhT_01210",
"RhT_00993",
"RhT_00329",
"RhT_00146",
"RhT_01062",
"RhT_01550",
"RhT_01672",
"RhT_00862",
"RhT_01638",
"RhT_01347",
"RhT_01315",
"RhT_00377",
"RhT_01004",
"RhT_01512",
"RhT_00502",
"RhT_01267",
"RhT_00643",
"RhT_00405",
"RhT_00189",
"RhT_01306",
"RhT_00319",
"RhT_00793",
"RhT_01308",
"RhT_00710",
"RhT_01397",
"RhT_00212",
"RhT_01034",
"RhT_01456",
"RhT_00176",
"RhT_01292",
"RhT_00786",
"RhT_01336",
"RhT_00795",
"RhT_01422",
"RhT_01252",
"RhT_00425",
"RhT_00494",
"RhT_00335",
"RhT_01420",
"RhT_00947",
"RhT_01297",
"RhT_00184",
"RhT_00493",
"RhT_00038",
"RhT_00315",
"RhT_01579",
"RhT_01421",
"RhT_01051",
"RhT_00369",
"RhT_01063",
"RhT_01581",
"RhT_00321",
"RhT_01502",
"RhT_01316",
"RhT_00873",
"RhT_01345",
"RhT_00330",
"RhT_00396",
"RhT_00632",
"RhT_00366",
"RhT_00948",
"RhT_01654",
"RhT_00586",
"RhT_01338",
"RhT_00889",
"RhT_01433",
"RhT_00486",
"RhT_00986",
"RhT_01426",
"RhT_00210",
"RhT_00732",
"RhT_01585",
"RhT_00705",
"RhT_00598",
"RhT_00765",
"RhT_00569",
"RhT_01625",
"RhT_00658",
"RhT_00574",
"RhT_00626",
"RhT_01666",
"RhT_00985",
"RhT_00987",
"RhT_00286",
"RhT_00971",
"RhT_01474",
"RhT_00429",
"RhT_01015",
"RhT_00186",
"RhT_00443",
"RhT_01310",
"RhT_00724",
"RhT_01435",
"RhT_01339",
"RhT_00714",
"RhT_01265",
"RhT_00633",
"RhT_00074",
"RhT_00290",
"RhT_01580",
"RhT_01380",
"RhT_00304",
"RhT_01496",
"RhT_01340",
"RhT_00223",
"RhT_00469",
"RhT_01261",
"RhT_00731",
"RhT_01152",
"RhT_00888",
"RhT_00417",
"RhT_00669",
"RhT_01359",
"RhT_01035",
"RhT_00949",
"RhT_00875",
"RhT_01570",
"RhT_01452",
"RhT_01521",
"RhT_00497",
"RhT_00155",
"RhT_01473",
"RhT_01296",
"RhT_01074",
"RhT_01285",
"RhT_00653",
"RhT_00941",
"RhT_01178",
"RhT_01262",
"RhT_00018",
"RhT_00154",
"RhT_01334",
"RhT_00563",
"RhT_01453",
"RhT_00980",
"RhT_00147",
"RhT_00361",
"RhT_01498",
"RhT_01414",
"RhT_00836",
"RhT_00650",
"RhT_01209",
"RhT_00421",
"RhT_00292",
"RhT_00302",
"RhT_00324",
"RhT_00921",
"RhT_01174",
"RhT_00410",
"RhT_00244",
"RhT_00037",
"RhT_00282",
"RhT_01419",
"RhT_01344",
"RhT_01165",
"RhT_01457",
"RhT_00604",
"RhT_01481",
"RhT_01413",
"RhT_01182",
"RhT_00984",
"RhT_01620",
"RhT_01288",
"RhT_01303",
"RhT_00524",
"RhT_01014",
"RhT_00222",
"RhT_01653",
"RhT_01075",
"RhT_01374",
"RhT_01642",
"RhT_01582",
"RhT_01439",
"RhT_01522",
"RhT_01222",
"RhT_00112",
"RhT_01104",
"RhT_01008",
"RhT_01107",
"RhT_01483",
"RhT_01111",
"RhT_01349",
"RhT_00991",
"RhT_00652",
"RhT_01476",
"RhT_01238",
"RhT_00635",
"RhT_00016",
"RhT_00697",
"RhT_00498",
"RhT_00291",
"RhT_01083",
"RhT_01415",
"RhT_00651",
"RhT_00200",
"RhT_00606",
"RhT_00877",
"RhT_01551",
"RhT_01430",
"RhT_00708",
"RhT_01651",
"RhT_00849",
"RhT_01279",
"RhT_01149",
"RhT_01478",
"RhT_01417",
"RhT_01412",
"RhT_00555",
"RhT_00853",
"RhT_00111",
"RhT_01377",
"RhT_01058",
"RhT_00994",
"RhT_01425",
"RhT_01097",
"RhT_00348",
"RhT_00882",
"RhT_00556",
"RhT_01423",
"RhT_00876",
"RhT_00202",
"RhT_00804",
"RhT_00686",
"RhT_01429",
"RhT_01305",
"RhT_01416",
"RhT_01379",
"RhT_00430",
"RhT_00002",
"RhT_01190",
"RhT_01432",
"RhT_00983",
"RhT_00211",
"RhT_00852",
"RhT_00975",
"RhT_00872",
"RhT_01431",
"RhT_01428",
"RhT_00201",
"RhT_00328",
"RhT_00871",
"RhT_00880",
"RhT_01514",
"RhT_01348",
"RhT_00054",
"RhT_01183",
"RhT_01084",
"RhT_00203"]
        #locus_highlight = conserved_chlamydiae
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
                                                                   draft_coordinates=False)





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

        circos_reference = gbk2circos.Circos_config(circos_files_reference["contigs"], chr_spacing_list)

        # add plus minus genes
        circos_reference.add_highlight(circos_files_reference["plus"], fill_color="grey_a1", r1="0.98r", r0="0.95r", href=href)
        circos_reference.add_highlight(circos_files_reference["minus"], fill_color="grey_a1", r1="0.95r", r0="0.92r", href=href)

        out_name = ''
        for accession in reference_accessions:
            out_name += accession

        out_file = "%s.svg" % out_name
        config_file = "%s.config" % out_name
        config_file_reference = os.path.join(out_directory, config_file)
        #out_file = os.path.join(out_directory, out_file)

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
        print cmd
        (stdout, stderr, return_code) = shell_command.shell_command(cmd)

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

        return r0

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
        return r1_minus, r0_minus

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
