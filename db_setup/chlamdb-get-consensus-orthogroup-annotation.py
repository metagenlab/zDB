#! /usr/bin/env python

def orthogroup_consensus_annotation(biodb, 
                                    COG=False, 
                                    interpro=False, 
                                    KO=False,
                                    gene=False):

    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    import operator

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select orthogroup_name,orthogroup_id from orthology_orthogroup;'
    group_name2group_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    if gene:
        print('Get orthogroup2gene')
        orthogroup2genes = biosql_own_sql_tables.orthogroup2gene(biodb)
        print('Get orthogroup2product')
        orthogroup2products = biosql_own_sql_tables.orthogroup2product(biodb)
    if COG:
        print('Get orthogroup2cog_id')
        orthogroup2cogs = biosql_own_sql_tables.orthogroup2cog_id(biodb)
    if interpro:
        print('Get orthogroup2pfam_id')
        orthogroup2pfam = biosql_own_sql_tables.orthogroup2pfam_id(biodb)
        print('Get orthogroup2interpro_id')
        orthogroup2interpro = biosql_own_sql_tables.orthogroup2interpro_id(biodb)
    if KO:
        print('Get orthogroup2ko_id')
        orthogroup2ko = biosql_own_sql_tables.orthogroup2ko_id(biodb)

    if gene:
        sql1 = 'create table if not exists orthology_orthogroup2gene (group_id INTEGER, `rank` INTEGER, count INTEGER, description TEXT);'
        sql2 = 'create table if not exists orthology_orthogroup2product (group_id INTEGER, `rank` INTEGER, count INTEGER, description TEXT);'
        server.adaptor.execute(sql1,)
        server.adaptor.execute(sql2,)
    if COG:
        sql3 = 'create table if not exists orthology_orthogroup2cog (group_id INTEGER, `rank` INTEGER, count INTEGER, COG_id INTEGER);'
        server.adaptor.execute(sql3,)
    if interpro:
        sql4 = 'create table if not exists orthology_orthogroup2pfam (group_id INTEGER, `rank` INTEGER, count INTEGER, signature_id INTEGER);'
        sql5 = 'create table if not exists orthology_orthogroup2interpro (group_id INTEGER, `rank` INTEGER, count INTEGER, interpro_id INTEGER);'
        server.adaptor.execute(sql4,)
        server.adaptor.execute(sql5,)
    if KO:
        sql6 = 'create table if not exists orthology_orthogroup2ko (group_id INTEGER, `rank` INTEGER, count INTEGER, ko_id INTEGER);'
        server.adaptor.execute(sql6,)

    template1 = 'insert into orthology_orthogroup2gene' + ' values (%s, %s, %s, %s)'
    template2 = 'insert into orthology_orthogroup2product' + ' values (%s, %s, %s, %s)'
    template3 = 'insert into orthology_orthogroup2cog' + ' values (%s, %s, %s, %s)'
    template4 = 'insert into orthology_orthogroup2pfam' + ' values (%s, %s, %s, %s)'
    template5 = 'insert into orthology_orthogroup2interpro' + ' values (%s, %s, %s, %s)'
    template6 = 'insert into orthology_orthogroup2ko' + ' values (%s, %s, %s, %s)'

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
        
        if COG:
            try:
                sorted_list = sorted(orthogroup2cogs[group].items(), key=operator.itemgetter(1), reverse=True)
                for n, row in enumerate(sorted_list):
                    cog, count = row
                    server.adaptor.execute(template3, [group_id, n, count, cog])
            except KeyError:
                continue
        if interpro:
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
        if KO:
            try:
                sorted_list = sorted(orthogroup2ko[group].items(), key=operator.itemgetter(1), reverse=True)
                for n, row in enumerate(sorted_list):
                    ko, count = row
                    server.adaptor.execute(template6, [group_id, n, count, ko])
            except KeyError:
                continue
    if gene:
        server.adaptor.commit()
        sql1 = 'create index ooggid on orthology_orthogroup2gene(group_id)'
        sql2 = 'create index oopgid on orthology_orthogroup2product(group_id)'
        server.adaptor.execute(sql1,)
        server.adaptor.execute(sql2,)
    if COG:
        sql3 = 'create index oocgid on orthology_orthogroup2cog(group_id)'
        sql7 = 'create index ooccid on orthology_orthogroup2cog(cog_id)'
        server.adaptor.execute(sql3,)
        server.adaptor.execute(sql7,)
    if interpro:
        sql4 = 'create index oopfgid on orthology_orthogroup2pfam(group_id)'
        sql5 = 'create index ooigid on orthology_orthogroup2interpro(group_id)'
        sql8 = 'create index oopsid on orthology_orthogroup2pfam(signature_id)'
        sql9 = 'create index ooipiid on orthology_orthogroup2interpro(interpro_id)'
        server.adaptor.execute(sql4,)
        server.adaptor.execute(sql5,)
        server.adaptor.execute(sql8,)
        server.adaptor.execute(sql9,)
    if KO:
        sql6 = 'create index ookgid on orthology_orthogroup2ko(group_id)'
        sql10 = 'create index ookkid on orthology_orthogroup2ko(ko_id)'  
        server.adaptor.execute(sql6,)
        server.adaptor.execute(sql10,)

    server.adaptor.commit()


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-g", '--gene', action="store_true", help="Gene and product")
    parser.add_argument("-c", '--COG', action="store_true", help="db name")
    parser.add_argument("-i", '--interpro', action="store_true", help="db name")
    parser.add_argument("-k", '--KO', action="store_true", help="db name")

    args = parser.parse_args()

    orthogroup_consensus_annotation(args.db_name,
                                    args.COG,
                                    args.interpro,
                                    args.KO,
                                    args.gene)
    
    # update config
    manipulate_biosqldb.update_config_table(args.db_name, "orthology_consensus_annotation")
