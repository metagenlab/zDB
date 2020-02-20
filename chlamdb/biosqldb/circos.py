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
                 locus_highlight2=[]):

        from chlamdb.biosqldb import manipulate_biosqldb
        from chlamdb.plots import gbk2circos
        import os
        from chlamdb.biosqldb import shell_command

        reference_accessions = []
        for reference_record in reference_records:
            #print "one record", reference_record
            reference_accessions.append(reference_record.id.split(".")[0])
        #print "reference_accession", reference_accessions


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

        klebsiella_missing = [
"group_4660",
"group_3913",
"group_2689",
"group_2688",
"group_3616",
"group_3324",
"group_3323",
"group_3322",
"group_136",
"group_2687",
"group_193",
"group_62",
"group_3610",
"group_3900",
"group_2691",
"group_3770",
"group_3940",
"group_3942",
"group_4414",
"group_3941",
"group_3301",
"group_3308",
"group_3302",
"group_3300",
"group_4472",
"group_3646",
"group_3628",
"group_2686",
"group_3295",
"group_3297",
"group_3296",
"group_3299",
"group_3298",
"group_100",
"group_3547",
"group_1339",
"group_2706",
"group_3838",
"group_3994",
"group_3995",
"group_4500",
"group_4042",
"group_4044",
"group_3938",
"group_3939",
"group_4260",
"group_4585",
"group_3914",
"group_3098",
"group_3619",
"group_3618",
"group_3617",
"group_83",
"group_88",
"group_3773",
"group_174",
"group_173",
"group_1345",
"group_60",
"group_4155",
"group_4154",
"group_4153",
"group_4038",
"group_3730",
"group_4713",
"group_4715",
"group_3101",
"group_3100",
"group_3103",
"group_3102",
"group_3105",
"group_3104",
"group_3107",
"group_3106",
"group_3109",
"group_3108",
"group_4594",
"group_2690",
"group_4342",
"group_3776",
"group_3771",
"group_4203",
"group_4405",
"group_3303",
"group_3307",
"group_3306",
"group_3305",
"group_237",
"group_3304",
"group_3624",
"group_3626",
"group_3627",
"group_3620",
"group_3621",
"group_3622",
"group_3623",
"group_3629",
"group_3607",
"group_4248",
"group_3551",
"group_3556",
"group_3112",
"group_4069",
"group_3099",
"group_2414",
"group_3110",
"group_3111",
"group_3726",
"group_4278",
"group_3772",
"group_4368",
"group_4343",
"group_3633",
"group_3791",
"group_3915",
"group_4286",
"group_3952",
"group_4732",
"group_4711",
"group_3625",
"group_3555",
"group_3839",
"group_4314",
"group_3853",
"group_3852",
"group_3851",
"group_3856",
"group_3855",
"group_3854",
"group_4397",
"group_3774",
"group_4474",
"group_4045",
"group_3879",
"group_3878",
"group_3912",
"group_3911",
"group_4795",
"group_4797",
"group_4401",
"group_4694",
"group_4793",
"group_4792",
"group_4794",
"group_4755",
"group_3947",
"group_4782",
"group_3997",


        ]











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
            sql = 'select locus_tag from blastnr.blastnr_%s t1 ' \
              ' inner join biosqldb.bioentry t2 on t1.query_bioentry_id=t2.bioentry_id ' \
              ' inner join biosqldb.biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
              ' inner join blastnr_blastnr_taxonomy t4 on t1.subject_taxid=t4.taxon_id ' \
              ' inner join custom_tables.locus2seqfeature_id_%s t5 ' \
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
