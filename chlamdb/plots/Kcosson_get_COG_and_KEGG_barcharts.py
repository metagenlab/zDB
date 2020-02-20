#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


locus_list = [

"KPN_RS08745",
"KPN_RS08800",
"KPN_RS08805",
"KPN_RS08810",
"KPN_RS08815",
"KPN_RS08820",
"KPN_RS08830",
"KPN_RS08840",
"KPN_RS08845",
"KPN_RS08850",
"KPN_RS08855",
"KPN_RS08860",
"KPN_RS09245",
"KPN_RS09250",
"KPN_RS09255",
"KPN_RS09830",
"KPN_RS09835",
"KPN_RS09840",
"KPN_RS09850",
"KPN_RS09855",
"KPN_RS09860",
"KPN_RS09865",
"KPN_RS09870",
"KPN_RS09875",
"KPN_RS09890",
"KPN_RS09895",
"KPN_RS09910",
"KPN_RS09925",
"KPN_RS09930",
"KPN_RS09935",
"KPN_RS09940",
"KPN_RS09945",
"KPN_RS09950",
"KPN_RS09955",
"KPN_RS09965",
"KPN_RS09970",
"KPN_RS09975",
"KPN_RS09980",
"KPN_RS09985",
"KPN_RS09990",
"KPN_RS09995",
"KPN_RS10000",
"KPN_RS10005",
"KPN_RS10010",
"KPN_RS10015",
"KPN_RS10020",
"KPN_RS10025",
"KPN_RS10030",
"KPN_RS10035",
"KPN_RS10040",
"KPN_RS10045",
"KPN_RS10050",
"KPN_RS10060",
"KPN_RS10065",
"KPN_RS10070",
"KPN_RS10075",
"KPN_RS10080",
"KPN_RS10085",
"KPN_RS10090",
"KPN_RS10095",
"KPN_RS10100",
"KPN_RS10105",
"KPN_RS10110",
"KPN_RS10115",
"KPN_RS10120",
"KPN_RS10130",
"KPN_RS10135",
"KPN_RS10140",
"KPN_RS10145",
"KPN_RS10150",
"KPN_RS10155",
"KPN_RS10160",
"KPN_RS10190",
"KPN_RS10195",
"KPN_RS10200",
"KPN_RS10205",
"KPN_RS10210",
"KPN_RS10215",
"KPN_RS10220",
"KPN_RS10225",
"KPN_RS10230",
"KPN_RS10235",
"KPN_RS10240",
"KPN_RS10765",
"KPN_RS12245",
"KPN_RS12250",
"KPN_RS12260",
"KPN_RS15910",
"KPN_RS15915",
"KPN_RS15920",
"KPN_RS15925",
"KPN_RS15930",
"KPN_RS15935",
"KPN_RS15940",
"KPN_RS15945",
"KPN_RS15950",
"KPN_RS15955",
"KPN_RS15960",
"KPN_RS15965",
"KPN_RS15970",
"KPN_RS15975",
"KPN_RS15980",
"KPN_RS15990",
"KPN_RS17990",
"KPN_RS18935",
"KPN_RS18965",
"KPN_RS23360",
"KPN_RS24085",
"KPN_RS24090",
"KPN_RS24095",
"KPN_RS24100",
"KPN_RS24105",
"KPN_RS24110",
"KPN_RS24115",
"KPN_RS24120",
"KPN_RS24125",
"KPN_RS24130",
"KPN_RS24135",
"KPN_RS24140",
"KPN_RS24145",
"KPN_RS24160",
"KPN_RS24165",
"KPN_RS24170",
"KPN_RS24175",
"KPN_RS24180",
"KPN_RS24185",
"KPN_RS24190",
"KPN_RS24195",
"KPN_RS24200",
"KPN_RS24205",
"KPN_RS24210",
"KPN_RS24215",
"KPN_RS24220",
"KPN_RS24225",
"KPN_RS24230",
"KPN_RS24235",
"KPN_RS24240",
"KPN_RS24245",
"KPN_RS24255",
"KPN_RS24260",
"KPN_RS24440",
"KPN_RS25290",
"KPN_RS25445",
"KPN_RS25450",
"KPN_RS25455"
]

from chlamdb.biosqldb import manipulate_biosqldb

server, db = manipulate_biosqldb.load_db("2017_03_30_kcosson")


def get_COG_data(locus_list):
    filter = '"'+'","'.join(locus_list)+'"'
    sql = 'select t4.code, t4.description, count(*) as n from COG.locus_tag2gi_hit_2017_03_30_kcosson t1 ' \
           'inner join COG_cog_names_2014 t2 on t1.COG_id=t2.COG_name ' \
           'inner join COG_cog_id2cog_category t3 on t2.COG_id=t3.COG_id ' \
           'inner join COG_code2category t4 on t3.category_id=t4.category_id' \
           ' where t1.locus_tag in (%s) group by t4.code,t4.description;' % filter
    data = server.adaptor.execute_and_fetchall(sql,)
    print sql
    for i in data:
        i = [str(n) for n in i]
        print '\t'.join(i)

def get_COG_data_ind(locus_list):
    data_list = []
    for locus in locus_list:
        sql = 'select locus_tag, COG_name, t2.description,t4.code, t4.description as n from COG.locus_tag2gi_hit_2017_03_30_kcosson t1 ' \
               'inner join COG_cog_names_2014 t2 on t1.COG_id=t2.COG_name ' \
               'inner join COG_cog_id2cog_category t3 on t2.COG_id=t3.COG_id ' \
               'inner join COG_code2category t4 on t3.category_id=t4.category_id' \
               ' where t1.locus_tag in ("%s");' % locus
        data = server.adaptor.execute_and_fetchall(sql,)
        #i = [n[0] for n in data]
        if len(data) != 0:
            for i in data:
                data_list.append(i)

    locus_tag2annot = {}
    for row in data_list:
        if row[0] not in locus_tag2annot:
            locus_tag2annot[row[0]] = [row[1], row[2], [row[3]], [row[4]]]
        else:
            locus_tag2annot[row[0]][2].append(row[3])
            locus_tag2annot[row[0]][3].append(row[4])
    #print locus_tag2annot
    for locus in locus_list:
        if locus not in locus_tag2annot:
            print
        else:
            #print locus_tag2annot[locus][2]
            print '%s\t%s\t%s\t%s\t%s' % (locus,
                                          locus_tag2annot[locus][0],
                                          locus_tag2annot[locus][1],
                                          ','.join(locus_tag2annot[locus][2]),
                                          ','.join(locus_tag2annot[locus][3]))

def get_KEGG_data(locus_list):
    filter = '"' + '","'.join(locus_list) + '"'
    sql = 'select t1.locus_tag,t1.ko_id, t2.definition, t4.* from enzyme.locus2ko_2017_03_30_kcosson t1 ' \
           'inner join enzyme_ko_annotation t2 on t1.ko_id=t2.ko_accession ' \
           'inner join enzyme_module2ko t3 on t2.ko_id=t3.ko_id ' \
           'inner join enzyme_kegg_module t4 on t3.module_id=t4.module_id ' \
           'where t1.locus_tag in (%s);' % filter
    sql = 'select module_sub_sub_cat,description, count(*) as n from enzyme.locus2ko_2017_03_30_kcosson t1 ' \
           'inner join enzyme_ko_annotation t2 on t1.ko_id=t2.ko_accession ' \
           'inner join enzyme_module2ko t3 on t2.ko_id=t3.ko_id ' \
           'inner join enzyme_kegg_module t4 on t3.module_id=t4.module_id ' \
           'where t1.locus_tag in (%s) group by module_sub_sub_cat, description;' % filter
    print sql
    data = server.adaptor.execute_and_fetchall(sql, )
    for i in data:
        i = [str(n) for n in i]
        print '\t'.join(i)

#get_COG_data(locus_list)
#get_KEGG_data(locus_list)
get_COG_data_ind(locus_list)