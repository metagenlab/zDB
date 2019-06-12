#!/usr/bin/env python

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import biosql_own_sql_tables

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", '--database_name', type=str, help="Database name")

    args = parser.parse_args()

    biosql_own_sql_tables.collect_genome_statistics(args.database_name)
