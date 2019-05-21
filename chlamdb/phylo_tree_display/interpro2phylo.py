#!/usr/bin/env python

from ete_motifs import interpro_tsv2pfam_data
from ete_motifs import draw_pfam_tree



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--tree',type=str,help="newick tree")
    parser.add_argument("-i", '--interpro_tsv', type=str, help="interpro tsv result table")
    parser.add_argument("-l", '--locus2organism', type=str, help="Tabulated file with: 'locustag\torganism_name' for leaves labels")

    args = parser.parse_args()

    locus2organism = {}
    with open(args.locus2organism, 'r') as f:
        for row in f:
            data = row.rstrip().split('\t')
            if len(data) != 2:
                raise IOError('expected a two column tsv for locus2organism table')
            else:
                locus2organism[data[0]] = data[1]
    interpro_data = interpro_tsv2pfam_data(args.interpro_tsv, locus2organism)

    t, ts, leaf_number = draw_pfam_tree(args.tree, interpro_data)

    path = 'pfam_tree.svg'

    if leaf_number < 10:
        leaf_number = 10
    ts.show_branch_support = True
    t.render(path, h=leaf_number*12, dpi=800, tree_style=ts)
