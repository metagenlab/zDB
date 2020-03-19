#!/usr/bin/env python

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import biosql_own_sql_tables
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", '--database_name', type=str, help="Database name")

    args = parser.parse_args()

    biosql_own_sql_tables.collect_genome_statistics(args.database_name)

    manipulate_biosqldb.update_config_table(args.database_name, "genome_statistics")