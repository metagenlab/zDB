#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-



def load_cog_names_table(table_file):
    '''
    create a table with distinct entries for cogs belonging to multiple categories (i.e NZ_CP007053  | COG2197 | TK)
    :param table_file:
    :return:
        '''
    import re
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db('chlamydia_04_16')
    sql = 'CREATE table COG_cog_names_2014 (COG_id varchar(100),' \
                                'functon varchar(10),' \
                                'name varchar(200), index COG_id(COG_id))'

    h = re.compile('^#.*')
    server.adaptor.execute_and_fetchall(sql,)

    with open(table_file, 'r') as f:
        for line in f:
            if re.match(h, line):
                continue
            else:
                data = line.rstrip().split('\t')
                cog = data[0]
                all_funct = list(data[1])
                name =  re.sub('"','', data[2])
                for one_category in all_funct:
                    sql = 'insert into COG_cog_names_2014 (COG_id, functon, name) values ("%s", "%s", "%s")' % (cog,
                                                                                                                 one_category,
                                                                                                                  name)
                    print sql
                    server.adaptor.execute(sql,)
                    server.commit()
load_cog_names_table('/home/trestan/stock/bio_databases/cog/cognames2003-2014.tab')