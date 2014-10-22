#! /usr/bin/env python


import sys
from Bio import GenBank
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from manipulate_biosqldb import load_db
from manipulate_biosqldb import query_yes_no

def import_gbk(gbk_name, server_name, db_name):
    parser = GenBank.FeatureParser()
    iterator = GenBank.Iterator(open(gbk_name), parser)
    records = [i for i in iterator]
    iterator = GenBank.Iterator(open(gbk_name), parser)
    if len(records)>1:
        print "genbank file contains more than one record (if the assembly contain more than one contig, contactenate contigs into a single record before importing it the the sql database"
        sys.exit()
    #With this iterator, the loading of the database is another one-liner:
    try:
        db_name.load(iterator)
    except:
        db_name.load(iterator)
        print gbk_name, "already into db?"
        pass
    server_name.adaptor.commit()



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    
    parser.add_argument("-g", '--gbk',type=str,help="gbk file", nargs='+')
    parser.add_argument("-d", '--db_name', type=str,help="db name")
    parser.add_argument("-r", '--remove_db', type=str,help="remove db")

    
    args = parser.parse_args()

    if not args.db_name:
        raise "db name?"
        
    server, db = load_db(args.db_name)

    if args.gbk:
        for gbk in args.gbk:
            print gbk
            import_gbk(gbk, server, db)
    
    if args.remove_db:
        if query_yes_no("Remove databse %s?" % args.remove_db):
            del server[args.remove_db]
            server.commit()

        
