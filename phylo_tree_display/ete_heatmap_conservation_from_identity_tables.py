#!/usr/bin/env python
'''

script to diplay results from get_chlam_genome_classif_seq.py => fasta => parwiseid_needle.py -m fasta => identity table
chose a reference and display percentages of identity with the other genomes

'''



import sys
from random import sample
from random import randint
from ete2 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace


import numpy as np
import colorsys
import matplotlib.colors as pcolors


def _get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        rgb = colorsys.hls_to_rgb(hue, lightness, saturation)

        print rgb
        colors.append(pcolors.rgb2hex(rgb))
    return colors




def organism2color(locus2data):

    organism_list = []
    for locus in locus2data:
        if locus2data[locus][0][2] not in organism_list:
            organism_list.append(locus2data[locus][0][2])
    colors = _get_colors(len(organism_list))

    return dict(zip(organism_list, colors))


def plot_heat_tree(gene2genome_id, tree_file, genes, accession2description=False, orgnames=False, taxid2st=False):
    '''
    Plot heatmap next to a tree. The order of the heatmap **MUST** be the same,
    as order of the leafs on the tree. The tree must be in the Newick format. If
    *output_file* is specified, then heat-tree will be rendered as a PNG,
    otherwise interactive browser will pop-up with your heat-tree.

    Parameters
    ----------
    heatmap_file: str
        Path to the heatmap file. The first row must have '#Names' as first
        element of the header.
            e.g. #Names, A, B, C, D
                row1, 2, 4, 0, 4
                row2, 4, 6, 2, -1

    tree_file: str
        Path to the tree file in Newick format. The leaf node labels should
        be the same as as row names in the heatmap file. E.g. row1, row2.

    output_file: str, optional
        If specified the heat-tree will be rendered in that file as a PNG image,
        otherwise interactive browser will pop-up. **N.B.** program will wait
        for you to exit the browser before continuing.
    '''
    import numpy

    from ete2.treeview.faces import add_face_to_node
    from ete2 import ClusterTree, TreeStyle, AttrFace, ProfileFace
    genome_list = []
    for gene in gene_list:
        for genome in gene2genome_id[gene]:
            if genome not in genome_list:
                genome_list.append(genome)


    print gene2genome_id['new_strain_Adk_identity']
    if orgnames:
        pass
        #taxid2organism = manipulate_biosqldb.taxon_id2genome_description(server, biodb, True)

    t1 = Tree(tree_file)
    #t.populate(8)
    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    # and set it as tree outgroup
    t1.set_outgroup(R)    # To operate with numbers efficiently

    family = ["16S", "23S"]

    genus = [
    "DnaA",
    "SucA",
    "protein_325",
    "FabI"
    ]
    species = [
    "RpoN",
    "FtsK",
    "PepF",
    "Adk",
    "HemL"]

    gene2cutoff = {
        '16S' : 92.5,
        '23S' : 91,
        'DnaA' : 70,
        'SucA' : 64,
        'protein_325' : 57,
        'FabI' : 78,
        "RpoN": 96,
        "FtsK": 98,
        "PepF": 96,
        "Adk": 95,
        "HemL": 95
    }

    file2category = {}
    file2gene = {}
    ordered_genes = family + genus + species
    ordered_files = []
    for gene in ordered_genes:
        for file in genes:
            if gene in file:
                ordered_files.append(file)
                file2gene[file] = gene
                if gene in genus:
                    file2category[file] = 'genus'
                elif gene in species:
                    file2category[file] = 'species'
                else:
                    file2category[file] = 'family'

    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=0.2, vmax=1) # map2count[map[0]][0]
    cmap_blue = cm.Blues
    cmap_red = cm.OrRd
    cmap_green = cm.BuGn

    m1 = cm.ScalarMappable(norm=norm, cmap=cmap_red)
    m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
    m3 = cm.ScalarMappable(norm=norm, cmap=cmap_green)

    leaf_number = 0
    for lf in t1.iter_leaves():
        leaf_number+=1

        lf.branch_vertical_margin = 0

        if orgnames:
            pass
        if taxid2st:
            try:
                st = taxid2st[lf.name]
            except:
                pass
        '''
        if accession2description:
            try:
                lf.name = accession2description[lf.name]
            except:
                pass
        '''

        for col, gene in enumerate(ordered_files):
            try:
                value = gene2genome_id[gene][lf.name]
            except:
                value = 0

            #print '#####', value, type(value)
            #print 'col', col, 'value', gene

            if leaf_number==1:
                n = TextFace('%s' % (file2gene[gene]), fsize=5)
                n.vt_align = 0
                n.hz_align = 1
                n.rotation= 270
                n.margin_top = 0
                n.margin_right = 0
                n.margin_left = 0
                n.margin_bottom = 0
                n.inner_background.color = "white"
                n.opacity = 1.
                lf.add_face(n, col, position="aligned")

            if value > 0:
                if len(str(value)) < 5:
                    n = TextFace(' %s   ' % value)
                else:
                    n = TextFace(' %s ' % value)
                n.margin_top = 0
                n.margin_right = 0
                n.margin_left = 0
                n.margin_bottom = 0
                if file2category[gene] == 'genus':
                    if value >= gene2cutoff[file2gene[gene]]:
                        n.inner_background.color = rgb2hex(m2.to_rgba(float(value/100))) #'#140718' #"#81BEF7"
                    else:
                        n.inner_background.color = 'grey'
                elif file2category[gene] == 'species':
                    if value >= gene2cutoff[file2gene[gene]]:
                        n.inner_background.color = rgb2hex(m1.to_rgba(float(value/100)))
                    else:
                        n.inner_background.color = 'grey'
                else:
                    if value >= gene2cutoff[file2gene[gene]]:
                        n.inner_background.color = rgb2hex(m3.to_rgba(float(value/100)))
                    else:
                        n.inner_background.color = 'grey'
                n.opacity = 1.
                lf.add_face(n, col, position="aligned")

            else:
                n = TextFace('   -   ')
                n.margin_top = 0
                n.margin_right = 0
                n.margin_left = 0
                n.margin_bottom = 0
                n.inner_background.color = "white"
                n.opacity = 1.
                lf.add_face(n, col, position="aligned")
    if accession2description:
        for lf in t1.iter_leaves():
            try:
                lf.name = accession2description[lf.name]
            except:
                pass
    #a = t1.render(output_file, dpi=800, h=i*15)
    return t1, leaf_number



def identity_tables2dico(reference_accession, table_list):
    '''
    for a given refernce genome, get the identity of each proteins in each of the other genomes

    :param reference_taxon:
    :param table_list: list of table with parwise identities
    :return: dictionnary gene2genome_id => with all identities
    '''

    import pandas as pd

    genome2gene_id = {}

    for table in table_list:
        gene = table.split('.')[0]
        data = pd.read_csv(table, sep="\t", index_col=0)

        try:
            genome2gene_id[gene] = dict(data[reference_accession])
        except:
            # gene not found
            genome2gene_id[gene] = {}

    return genome2gene_id

if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--tree',
                        type=str,
                        help="newick tree")
    parser.add_argument("-m",'--identity_tables',
                        nargs='+',
                        help="identity tables (tab file with gene names as columns and genomes as rows)")
    parser.add_argument("-s",'--mlst',
                        type=str,
                        help="mlst file",
                        default=False)
    parser.add_argument("-a",'--accession2description',
                        default=False,
                        type=str,
                        help="tab file with accessions and corresponding descriptions for leaf labels")
    parser.add_argument("-g",'--gbk_files',
                        default=False,
                        nargs='+',
                        help="genbank files to get accession2description conversion")
    parser.add_argument("-r",'--reference',
                        type=str,
                        help="reference_accession")

    args = parser.parse_args()

    print len(args.identity_tables)

    gene2genome_id = identity_tables2dico(args.reference, args.identity_tables)
    gene_list = gene2genome_id.keys()

    genome_list = []
    for gene in gene_list:
        for genome in gene2genome_id[gene]:
            if genome not in genome_list:
                genome_list.append(genome)
    print genome_list
    print gene_list

    if args.accession2description:
        accession2description = {}
        with open(args.accession2description, 'r') as f:
            for row in f:
                data = row.rstrip().split('\t')
                accession2description[data[0]] = data[1]
    elif args.gbk_files:
        import gbk2accessiontodefinition
        accession2description = gbk2accessiontodefinition.get_coressp(args.gbk_files)
    else:
        accession2description = False

    if args.mlst:
        taxid2st = {}
        with open(args.mlst) as f:
            for n, row in enumerate(f):
                if n==0:
                    continue
                else:
                    data = row.rstrip().split('\t')
                    taxid2st[data[0]] = data[2]
    else:
        taxid2st = False

    t, n = plot_heat_tree(gene2genome_id, args.tree, gene_list, accession2description, taxid2st)
    t.render("test.svg")