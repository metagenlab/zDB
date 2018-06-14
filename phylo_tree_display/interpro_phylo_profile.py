#!/usr/bin/env python

def parse_interpro_data(interpro_file, analysis="Pfam"):


def interpro2ete_tree(ete2_tree, interproscan_table):




if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--table',type=str,help="Interproscan table")
    parser.add_argument("-p", '--phylogeny', help="Phylogenetic tree")
    parser.add_argument("-a", '--analysis', default='Pfam', help="Interproscan analysis to display")

    args = parser.parse_args()
    #plot_phylum_counts()


