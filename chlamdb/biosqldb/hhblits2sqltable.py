#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def parse_a3m(a3m_file):
    import Bio.SeqIO
    import A3MIO
    with open(a3m_file, 'r') as f:
        for n,line in enumerate(f):
            try:
                data = line[36:].rstrip().split() + [line[0:35].split()[1].split('|')[1]]
                print data, len(data)
                if n==20:
                    break
            except:
                pass

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import sys
    import json

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_blast', type=str, help="input blast tab files", nargs='+')
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-t", '--create_tables', action='store_true', help="Create SQL tables")
    parser.add_argument("-p", '--n_procs', type=int, help="Number of threads to use (default=8)", default=8)
    parser.add_argument("-c", '--clean_tables', action='store_true', help="delete all sql tables")

    parse_a3m("/home/trestan/work/dev/genomic_dbs/chlamydia_04_16b/test_hhblits/AMV16246_1.hhr")

