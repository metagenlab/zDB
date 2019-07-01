#! /usr/bin/env python

def orthogroup_consensus_annotation(biodb):

    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    import operator

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select orthogroup_name,orthogroup_id from orthology.orthogroup_%s;' % biodb
    group_name2group_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    print('Get orthogroup2gene')
    orthogroup2genes = biosql_own_sql_tables.orthogroup2gene(biodb)
    print('Get orthogroup2product')
    orthogroup2products = biosql_own_sql_tables.orthogroup2product(biodb)
    print('Get orthogroup2cog_id')
    orthogroup2cogs = biosql_own_sql_tables.orthogroup2cog_id(biodb)
    print('Get orthogroup2pfam_id')
    orthogroup2pfam = biosql_own_sql_tables.orthogroup2pfam_id(biodb)
    print('Get orthogroup2interpro_id')
    orthogroup2interpro = biosql_own_sql_tables.orthogroup2interpro_id(biodb)
    print('Get orthogroup2ko_id')
    orthogroup2ko = biosql_own_sql_tables.orthogroup2ko_id(biodb)

    sql1 = 'create table if not exists orthology.orthogroup2gene_%s (group_id INTEGER, rank INTEGER, count INTEGER, description TEXT);' % biodb
    sql2 = 'create table if not exists orthology.orthogroup2product_%s (group_id INTEGER, rank INTEGER, count INTEGER, description TEXT);' % biodb
    sql3 = 'create table if not exists orthology.orthogroup2cog_%s (group_id INTEGER, rank INTEGER, count INTEGER, COG_id INTEGER);' % biodb
    sql4 = 'create table if not exists orthology.orthogroup2pfam_%s (group_id INTEGER, rank INTEGER, count INTEGER, signature_id INTEGER);' % biodb
    sql5 = 'create table if not exists orthology.orthogroup2interpro_%s (group_id INTEGER, rank INTEGER, count INTEGER, interpro_id INTEGER);' % biodb
    sql6 = 'create table if not exists orthology.orthogroup2ko_%s (group_id INTEGER, rank INTEGER, count INTEGER, ko_id INTEGER);' % biodb

    server.adaptor.execute(sql1,)
    server.adaptor.execute(sql2,)
    server.adaptor.execute(sql3,)
    server.adaptor.execute(sql4,)
    server.adaptor.execute(sql5,)
    server.adaptor.execute(sql6,)

    template1 = 'insert into orthology.orthogroup2gene_%s' % biodb + ' values (%s, %s, %s, %s)'
    template2 = 'insert into orthology.orthogroup2product_%s' % biodb + ' values (%s, %s, %s, %s)'
    template3 = 'insert into orthology.orthogroup2cog_%s' % biodb + ' values (%s, %s, %s, %s)'
    template4 = 'insert into orthology.orthogroup2pfam_%s' % biodb + ' values (%s, %s, %s, %s)'
    template5 = 'insert into orthology.orthogroup2interpro_%s' % biodb + ' values (%s, %s, %s, %s)'
    template6 = 'insert into orthology.orthogroup2ko_%s' % biodb + ' values (%s, %s, %s, %s)'

    for n, group in enumerate(group_name2group_id):
        print(n, group)
        group_id = group_name2group_id[group]
        try:
            sorted_list = sorted(orthogroup2genes[group].items(), key=operator.itemgetter(1), reverse=True)
            for n, row in enumerate(sorted_list):
                gene, count = row
                server.adaptor.execute(template1, [group_id, n, count, gene])
        except KeyError:
            continue
        try:
            sorted_list = sorted(orthogroup2products[group].items(), key=operator.itemgetter(1), reverse=True)
            for n, row in enumerate(sorted_list):
                product, count = row
                server.adaptor.execute(template2, [group_id, n, count, product])
        except KeyError:
            continue
        try:
            sorted_list = sorted(orthogroup2cogs[group].items(), key=operator.itemgetter(1), reverse=True)
            for n, row in enumerate(sorted_list):
                cog, count = row
                server.adaptor.execute(template3, [group_id, n, count, cog])
        except KeyError:
            continue
        try:
            sorted_list = sorted(orthogroup2pfam[group].items(), key=operator.itemgetter(1), reverse=True)
            for n, row in enumerate(sorted_list):
                pfam, count = row
                server.adaptor.execute(template4, [group_id, n, count, pfam])
        except KeyError:
            continue
        try:
            sorted_list = sorted(orthogroup2interpro[group].items(), key=operator.itemgetter(1), reverse=True)
            for n, row in enumerate(sorted_list):
                interpro, count = row
                server.adaptor.execute(template5, [group_id, n, count, interpro])
        except KeyError:
            continue
        try:
            sorted_list = sorted(orthogroup2ko[group].items(), key=operator.itemgetter(1), reverse=True)
            for n, row in enumerate(sorted_list):
                ko, count = row
                server.adaptor.execute(template6, [group_id, n, count, ko])
        except KeyError:
            continue

    server.adaptor.commit()
    sql1 = 'create index gg on orthology.orthogroup2gene_%s(group_id)' % biodb
    sql2 = 'create index gpr on orthology.orthogroup2product_%s(group_id)' % biodb
    sql3 = 'create index gc on orthology.orthogroup2cog_%s(group_id)' % biodb
    sql4 = 'create index gpf on orthology.orthogroup2pfam_%s(group_id)' % biodb
    sql5 = 'create index gpf on orthology.orthogroup2interpro_%s(group_id)' % biodb
    sql6 = 'create index gpf on orthology.orthogroup2ko_%s(group_id)' % biodb
    sql7 = 'create index id1 on orthology.orthogroup2cog_%s(cog_id)' % biodb
    sql8 = 'create index id2 on orthology.orthogroup2pfam_%s(signature_id)' % biodb
    sql9 = 'create index id3 on orthology.orthogroup2interpro_%s(interpro_id)' % biodb
    sql10 = 'create index id4 on orthology.orthogroup2ko_%s(ko_id)' % biodb
    server.adaptor.execute(sql1,)
    server.adaptor.execute(sql2,)
    server.adaptor.execute(sql3,)
    server.adaptor.execute(sql4,)
    server.adaptor.execute(sql5,)
    server.adaptor.execute(sql6,)
    server.adaptor.execute(sql7,)
    server.adaptor.execute(sql8,)
    server.adaptor.execute(sql9,)
    server.adaptor.execute(sql10,)

    server.adaptor.commit()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    orthogroup_consensus_annotation(args.db_name)
