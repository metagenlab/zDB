#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def biodb2pairwise_dist_phylogenies(biodb):

    from ete2 import Tree
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select phylogeny from phylogenies' % biodb
    all_phylogenies = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    sql = 'select seqfeature_id, taxon_id from custom_tables_locus2seqfeature_id' % biodb
    seqfeature_id2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'create table if not exists phylogenies.pairwise_leaf_dist_%s (taxon_1 INT, ' \
          ' taxon_2 INT, ' \
          ' seqfeature_1 INT, ' \
          ' seqfeature_2 INT, ' \
          ' distance FLOAT, ' \
          ' index taxon_1(taxon_1), ' \
          ' index taxon_2(taxon_2), ' \
          ' index seqfeature_1(seqfeature_1),' \
          ' index seqfeature_2(seqfeature_2), ' \
          ' index distance(distance))' % (biodb)

    server.adaptor.execute(sql,)
    for c, phylogeny in enumerate(all_phylogenies):

        t = Tree(phylogeny)
        leaf_list = [i for i in t.iter_leaves()]
        print '%s / %s -- %s' % (c, len(all_phylogenies), len(leaf_list))
        for n, l1 in enumerate(leaf_list):
            #print "%s / %s" % (n, len(leaf_list))
            for l2 in leaf_list[n+1:]:
                seqid_1 = locus_tag2seqfeature_id[str(l1.name)]
                seqid_2 = locus_tag2seqfeature_id[str(l2.name)]
                taxon_1 = seqfeature_id2taxon_id[str(seqid_1)]
                taxon_2 = seqfeature_id2taxon_id[str(seqid_2)]

                sql = 'insert into phylogenies.pairwise_leaf_dist_%s values (%s, %s, %s, %s, %s)' % (biodb,
                                                                                                              taxon_1,
                                                                                                              taxon_2,
                                                                                                              seqid_1,
                                                                                                              seqid_2,
                                                                                                              l1.get_distance(l2))
                server.adaptor.execute(sql,)
        server.commit()


def get_normalized_pairwise_dist(biodb, sub=False):

    from chlamdb.biosqldb import manipulate_biosqldb
    import numpy

    server, db = manipulate_biosqldb.load_db(biodb)
    sql = 'select distinct taxon_id from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t2.name="%s"' % biodb

    taxon_id_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    if sub:
        sql = 'create table if not exists phylogenies.pairwise_leaf_dist_norm_sub_%s (taxon_1 INT, ' \
              ' taxon_2 INT, ' \
              ' seqfeature_1 INT, ' \
              ' seqfeature_2 INT, ' \
              ' distance_median_norm FLOAT, ' \
              ' index taxon_1(taxon_1), ' \
              ' index taxon_2(taxon_2), ' \
              ' index seqfeature_1(seqfeature_1),' \
              ' index seqfeature_2(seqfeature_2), ' \
              ' index distance_median_norm(distance_median_norm))' % (biodb)
    else:
        sql = 'create table if not exists phylogenies.pairwise_leaf_dist_norm_%s (taxon_1 INT, ' \
              ' taxon_2 INT, ' \
              ' seqfeature_1 INT, ' \
              ' seqfeature_2 INT, ' \
              ' distance_median_norm FLOAT, ' \
              ' index taxon_1(taxon_1), ' \
              ' index taxon_2(taxon_2), ' \
              ' index seqfeature_1(seqfeature_1),' \
              ' index seqfeature_2(seqfeature_2), ' \
              ' index distance_median_norm(distance_median_norm))' % (biodb)

    server.adaptor.execute(sql,)

    for n, taxon_id_1 in enumerate(taxon_id_list):
        for taxon_id_2 in taxon_id_list[n+1:]:
            print 'taxons: %s-%s' % (taxon_id_1, taxon_id_2)

            sql = 'select taxon_1, taxon_2, seqfeature_1, seqfeature_2, distance from phylogenies.pairwise_leaf_dist_%s ' \
                  ' where taxon_1=%s and taxon_2=%s UNION select taxon_2, taxon_1, seqfeature_2, seqfeature_1, distance ' \
                  ' from phylogenies.pairwise_leaf_dist_%s  where taxon_1=%s and taxon_2=%s' % (biodb,
                                                                   taxon_id_1,
                                                                   taxon_id_2,
                                                                   biodb,
                                                                   taxon_id_2,
                                                                   taxon_id_1)

            data = server.adaptor.execute_and_fetchall(sql,)
            distance_list = [float(i[4]) for i in data]
            dist_median = numpy.median(distance_list)
            dist_mean = numpy.mean(distance_list)
            print dist_mean, dist_median
            for row in data:
                taxon_1 = row[0]
                taxon_2 = row[1]
                seqfeature_1 = row[2]
                seqfeature_2 = row[3]
                if sub:
                    norm_dist = float(row[4])-dist_median
                    sql = 'insert into phylogenies.pairwise_leaf_dist_norm_sub_%s values (%s, %s, %s, %s, %s)' % (biodb,
                                                                                                                taxon_1,
                                                                                                                taxon_2,
                                                                                                                seqfeature_1,
                                                                                                                seqfeature_2,
                                                                                                                norm_dist)
                else:
                    norm_dist = float(row[4])/dist_median
                    sql = 'insert into phylogenies.pairwise_leaf_dist_norm_%s values (%s, %s, %s, %s, %s)' % (biodb,
                                                                                                                taxon_1,
                                                                                                                taxon_2,
                                                                                                                seqfeature_1,
                                                                                                                seqfeature_2,
                                                                                                                norm_dist)
                server.adaptor.execute(sql,)
        server.commit()




def get_orthogroup_median_dist(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    import numpy
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table if not exists phylogenies.orthogroup2median_distance_%s (orthogroup varchar(200),' \
          ' median_norm_dist FLOAT, index orthogroup(orthogroup))' % biodb
    print sql
    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select orthogroup, locus_tag from orthology_detail' % biodb

    data = server.adaptor.execute_and_fetchall(sql,)
    orthogroup2locus_list = {}
    for row in data:
        if row[0] not in orthogroup2locus_list:
            orthogroup2locus_list[row[0]] = [row[1]]
        else:
            orthogroup2locus_list[row[0]].append(row[1])

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb

    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for n, one_group in enumerate(orthogroup2locus_list):
        print '%s - %s / %s' % (one_group, n, len(orthogroup2locus_list))

        locus_list = orthogroup2locus_list[one_group]

        if len(locus_list) < 4:
            continue

        dist_list = []
        for l1, locus_1 in enumerate(locus_list):
            for l2, locus_2 in enumerate(locus_list[l1+1:]):

                sql = 'select distance_median_norm from phylogenies.pairwise_leaf_dist_norm_sub_%s ' \
                      ' where (seqfeature_1 =%s and seqfeature_2=%s or ' \
                      '  seqfeature_1 =%s and seqfeature_2=%s)' % (biodb,
                                                                       locus_tag2seqfeature_id[locus_1],
                                                                       locus_tag2seqfeature_id[locus_2],
                                                                   locus_tag2seqfeature_id[locus_2],
                                                                   locus_tag2seqfeature_id[locus_1]
                                                                   )
                #print sql

                try:
                    #print one_group, server.adaptor.execute_and_fetchall(sql,)[0][0]
                    dist_list.append(float(server.adaptor.execute_and_fetchall(sql,)[0][0]))
                except:
                    pass
                    #print 'problem!', sql
        #print 'list', dist_list

        if len(dist_list) > 0:
            median_dist = numpy.median(dist_list)
            sql = 'insert into phylogenies.orthogroup2median_distance_%s values ("%s", %s)' % (biodb,
                                                                                                        one_group,
                                                                                                        str(median_dist))
            print sql
            server.adaptor.execute(sql,)
            server.commit()

#biodb2pairwise_dist_phylogenies('chlamydia_04_16')
#get_normalized_pairwise_dist('chlamydia_04_16', sub=True)
#get_orthogroup_median_dist('chlamydia_04_16')

def get_most_conserved(biodb):

    grp_list = ["group_619",
"group_613",
"group_429",
"group_428",
"group_427",
"group_422",
"group_421",
"group_328",
"group_324",
"group_325",
"group_326",
"group_327",
"group_320",
"group_321",
"group_322",
"group_323",
"group_128",
"group_129",
"group_125",
"group_122",
"group_123",
"group_121",
"group_695",
"group_696",
"group_711",
"group_710",
"group_715",
"group_717",
"group_716",
"group_95",
"group_94",
"group_90",
"group_232",
"group_233",
"group_834",
"group_830",
"group_831",
"group_476",
"group_477",
"group_471",
"group_472",
"group_473",
"group_500",
"group_505",
"group_506",
"group_507",
"group_508",
"group_509",
"group_286",
"group_285",
"group_283",
"group_975",
"group_1026",
"group_1024",
"group_1023",
"group_231",
"group_218",
"group_581",
"group_621",
"group_622",
"group_299",
"group_290",
"group_291",
"group_293",
"group_297",
"group_68",
"group_66",
"group_65",
"group_62",
"group_60",
"group_61",
"group_315",
"group_317",
"group_316",
"group_311",
"group_313",
"group_312",
"group_319",
"group_318",
"group_1113",
"group_1110",
"group_157",
"group_156",
"group_154",
"group_152",
"group_151",
"group_150",
"group_479",
"group_807",
"group_806",
"group_808",
"group_399",
"group_398",
"group_395",
"group_396",
"group_391",
"group_392",
"group_531",
"group_532",
"group_535",
"group_537",
"group_275",
"group_271",
"group_1017",
"group_1015",
"group_749",
"group_744",
"group_742",
"group_740",
"group_882",
"group_228",
"group_920",
"group_674",
"group_348",
"group_349",
"group_342",
"group_343",
"group_340",
"group_346",
"group_347",
"group_344",
"group_345",
"group_1123",
"group_1126",
"group_1125",
"group_1128",
"group_1129",
"group_188",
"group_189",
"group_184",
"group_180",
"group_182",
"group_183",
"group_397",
"group_39",
"group_34",
"group_36",
"group_31",
"group_30",
"group_33",
"group_32",
"group_445",
"group_447",
"group_446",
"group_440",
"group_443",
"group_449",
"group_859",
"group_460",
"group_104",
"group_105",
"group_106",
"group_100",
"group_101",
"group_109",
"group_221",
"group_223",
"group_227",
"group_738",
"group_565",
"group_560",
"group_957",
"group_279",
"group_601",
"group_604",
"group_418",
"group_419",
"group_412",
"group_413",
"group_410",
"group_411",
"group_416",
"group_414",
"group_415",
"group_96",
"group_339",
"group_338",
"group_333",
"group_332",
"group_331",
"group_337",
"group_336",
"group_335",
"group_334",
"group_138",
"group_134",
"group_137",
"group_136",
"group_133",
"group_132",
"group_686",
"group_684",
"group_688",
"group_267",
"group_266",
"group_262",
"group_88",
"group_89",
"group_82",
"group_493",
"group_495",
"group_828",
"group_519",
"group_518",
"group_516",
"group_515",
"group_514",
"group_510",
"group_1030",
"group_1031",
"group_1032",
"group_1033",
"group_1036",
"group_1038",
"group_1068",
"group_1062",
"group_763",
"group_207",
"group_206",
"group_205",
"group_203",
"group_202",
"group_201",
"group_597",
"group_593",
"group_592",
"group_599",
"group_909",
"group_904",
"group_906",
"group_652",
"group_654",
"group_970",
"group_973",
"group_974",
"group_280",
"group_58",
"group_52",
"group_57",
"group_56",
"group_55",
"group_361",
"group_362",
"group_363",
"group_364",
"group_365",
"group_366",
"group_367",
"group_368",
"group_161",
"group_165",
"group_468",
"group_462",
"group_461",
"group_467",
"group_466",
"group_465",
"group_464",
"group_1029",
"group_1069",
"group_1061",
"group_1060",
"group_1065",
"group_548",
"group_547",
"group_934",
"group_932",
"group_939",
"group_669",
"group_668",
"group_1124",
"group_258",
"group_254",
"group_257",
"group_250",
"group_253",
"group_494",
"group_358",
"group_351",
"group_350",
"group_357",
"group_356",
"group_193",
"group_191",
"group_190",
"group_197",
"group_196",
"group_211",
"group_212",
"group_22",
"group_23",
"group_24",
"group_430",
"group_431",
"group_432",
"group_434",
"group_435",
"group_436",
"group_437",
"group_438",
"group_439",
"group_843",
"group_246",
"group_113",
"group_110",
"group_117",
"group_116",
"group_119",
"group_1074",
"group_1078",
"group_1079",
"group_706",
"group_704",
"group_705",
"group_490",
"group_496",
"group_497",
"group_448",
"group_576",
"group_573",
"group_969",
"group_962",
"group_269",
"group_268",
"group_444",
"group_442",
"group_637",
"group_74",
"group_77",
"group_409",
"group_408",
"group_401",
"group_400",
"group_403",
"group_402",
"group_405",
"group_404",
"group_407",
"group_406",
"group_394",
"group_307",
"group_305",
"group_301",
"group_308",
"group_309",
"group_393",
"group_142",
"group_145",
"group_147",
"group_1086",
"group_1081",
"group_1080",
"group_1083",
"group_1088",
"group_489",
"group_488",
"group_810",
"group_558",
"group_584",
"group_389",
"group_529",
"group_526",
"group_1000",
"group_773",
"group_770",
"group_776",
"group_893",
"group_239",
"group_481",
"group_483",
"group_482",
"group_485",
"group_484",
"group_487",
"group_486",
"group_919",
"group_645",
"group_643",
"group_642",
"group_815",
"group_49",
"group_42",
"group_43",
"group_44",
"group_45",
"group_993",
"group_991",
"group_990",
"group_997",
"group_996",
"group_994",
"group_379",
"group_378",
"group_377",
"group_376",
"group_375",
"group_374",
"group_373",
"group_372",
"group_371",
"group_370",
"group_1130",
"group_171",
"group_170",
"group_174",
"group_177",
"group_456",
"group_457",
"group_454",
"group_455",
"group_452",
"group_453",
"group_450",
"group_451",
"group_458",
"group_459",
"group_867",
"group_865",
"group_860",
"group_384",
"group_385",
"group_382",
"group_383",
"group_380",
"group_381",
"group_222",
"group_1076",
"group_1077",
"group_1072",
"group_1073",
"group_727",
"group_550",
"group_557",
"group_788",
"group_783",
"group_2",
"group_1",
"group_6",
"group_5",
"group_4",
"group_940",
"group_941",
"group_242",
"group_241",
"group_248"]
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    filter = '"'+'","'.join(grp_list)+'"'

    sql = 'select A.*, B.locus_tag, B.gene, B.product from (select * from phylogenies.orthogroup2median_distance_%s where orthogroup in (%s) ' \
          ' order by median_norm_dist limit 550) A inner join orthology_detail B on A.orthogroup=B.orthogroup ' \
          ' where locus_tag like "%%%%WCW%%%%" order by median_norm_dist' % (biodb, filter, biodb)
    data = server.adaptor.execute_and_fetchall(sql)

    for row in data:
        #print row
        row = [str(i) for i in row]
        print '\t'.join(row)

get_most_conserved('chlamydia_04_16')