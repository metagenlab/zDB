#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys
from random import sample
from random import randint
from ete2 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace, faces, NodeStyle


import numpy as np
import colorsys
import matplotlib.colors as pcolors
import matplotlib.cm as cm
from matplotlib.colors import rgb2hex
import matplotlib as mpl

norm = mpl.colors.Normalize(vmin=80, vmax=100)
cmap_red = cm.OrRd
cmap_blue = cm.Blues
cmap_green = cm.Greens
m_red = cm.ScalarMappable(norm=norm, cmap=cmap_red)
m_blue = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
m_green = cm.ScalarMappable(norm=norm, cmap=cmap_green)



def parse_mlst_results(mlst_file):
    '''
    Parse mlst result from script wrote by torsen seeman
    return accession2ST_type
    :return:
    '''

    accession2st_type = {}
    with open(mlst_file, 'r') as f:
        for line in f:
            line = line.rstrip().split('\t')
            accession2st_type[line[0]] = line[2]
    return accession2st_type


def remove_blast_redundancy(blast_file_list, check_overlap=True):

    blast2data = {}
    queries = []

    # keep only best hit
    # todo make it proper
    for one_blast_file in blast_file_list:

        result_handle = open(one_blast_file, 'r')
        best_hit_handle = open("niq_best.tmp", 'w')

        hit_list = []
        for line in result_handle:
            if line.split('\t')[0] in hit_list:
                continue
            else:
                hit_list.append(line.split('\t')[0])
                best_hit_handle.write(line)
        best_hit_handle.close()


        with open("niq_best.tmp", 'r') as f:
            for line in f:
                line = line.split('\t')
                if line[1] not in blast2data:
                    blast2data[line[1]] = {}
                    blast2data[line[1]][line[0]] = [float(line[2]), int(line[8]), int(line[9])]
                else:
                     blast2data[line[1]][line[0]] = [float(line[2]), int(line[8]), int(line[9])]
                if line[0] not in queries:
                    queries.append(line[0])
    #print 'raw dico'
    #print blast2data

    print 'preparing blast dico'
    if check_overlap:
        for one_blast in blast2data.keys():
            for ref_gene in blast2data[one_blast].keys():

                for query_gene in blast2data[one_blast].keys():
                    overlap = False
                    if ref_gene == query_gene:
                        continue
                    # check if position is overlapping
                    try:
                        sorted_coordinates = sorted(blast2data[one_blast][ref_gene][1:3])
                        if blast2data[one_blast][query_gene][1] <= sorted_coordinates[1] and blast2data[one_blast][query_gene][1]>= sorted_coordinates[0]:
                            print 'Overlaping locations!'
                            print one_blast, ref_gene, blast2data[one_blast][ref_gene]
                            print one_blast, query_gene, blast2data[one_blast][query_gene]
                            overlap =True
                        sorted_coordinates = sorted(blast2data[one_blast][query_gene][1:3])
                        if blast2data[one_blast][ref_gene][1] <= sorted_coordinates[1] and blast2data[one_blast][ref_gene][1]>= sorted_coordinates[0]:
                            print 'Overlapping locations!'
                            print one_blast, ref_gene, blast2data[one_blast][ref_gene]
                            print one_blast, query_gene, blast2data[one_blast][query_gene]
                            overlap =True
                        if overlap:
                            if float(blast2data[one_blast][ref_gene][0]) > float(blast2data[one_blast][query_gene][0]):
                                del blast2data[one_blast][query_gene]
                                print 'removing', query_gene
                            else:
                                del blast2data[one_blast][ref_gene]
                                print 'removing', ref_gene
                                break
                    except KeyError:
                        pass
                        #print 'colocation already resolved:', one_blast, query_gene, ref_gene
    return blast2data, queries



def plot_blast_result(tree_file,
                      blast_result_file_list,
                      id2description,
                      id2mlst,check_overlap,
                      ordered_queries,
                      id_cutoff=80,
                      reference_accession='-',
                      accession2hit_filter=False,
                      show_identity_values=True):
    '''
    Projet Staph aureus PVL avec Laure Jaton
    Script pour afficher une phylogénie et la conservation de facteurs de virulence côte à côte
    Nécessite résultats MLST, ensemble des résultats tblastn (facteurs de virulence vs chromosomes),
    ainsi qu'une correspondance entre les accession des génomes et les noms qu'on veut afficher dans la phylogénie. Icemn
    pour les identifiants molis des patients, on les remplace par CHUV n.
    :param tree_file: phylogénie au format newick avec identifiants correspondants à tous les dico utilisés
    :param blast_result_file_list: résultats tblastn virulence factors vs chromosome (seulement best blast)
    :param id2description: identifiants génome utiisé dans l'arbre et description correspondante (i.e S aureus Newman)
    :param id2mlst: identitifiants arbre 2 S. aureus ST type
    :return:
    '''

    blast2data, queries = remove_blast_redundancy(blast_result_file_list, check_overlap)

    print 'blast dico ok, delete columns with no matches'
    queries_count = {}

    for query in queries:
        queries_count[query] = 0
        for one_blast in blast2data:
            if query in blast2data[one_blast]:
                #print blast2data[one_blast][query]
                if float(blast2data[one_blast][query][0]) > id_cutoff:
                    queries_count[query]+=1
                else:
                    del blast2data[one_blast][query]
    print 'query counts', len(queries_count)
    print queries_count
    for query in queries:
        print "%s\t%s" % (query, queries_count[query])
        if queries_count[query] == 0:
            queries.pop(queries.index(query))


    

    print 'delete columns with no matches ok'

    '''             
    rm_genes = ['selv','spsmA1','psmB1','psmB2','ses','set','sel','selX','sek','sel2','LukF', 'LukM', 'hly', 'hld'
        , 'hlgA', 'hlgB', 'hlgC', 'sed', 'sej', 'ser', 'selq1', 'sec3', 'sek2', 'seq2', 'lukD', 'lukE']
    #rm_genes = ['icaR','icaA','icaB','icaC','icaD', 'sdrF', 'sdrH']

    for gene in rm_genes:
        queries.pop(queries.index(gene))
    '''
    #queries = ['selv']
    t1 = Tree(tree_file)
    tss = TreeStyle()
    #tss.show_branch_support = True
    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()

    ordered_queries_filtered = []
    for query in ordered_queries:
        hit_count = 0
        for lf2 in t1.iter_leaves():
            try:
                tmpidentity = blast2data[lf2.name][query][0]
                if float(tmpidentity)>float(id_cutoff):
                    hit_count+=1
            except:
                continue
        if hit_count>0:
            ordered_queries_filtered.append(query)
            #print 'skippink-----------'
    
    head = True
    print 'drawing tree'
    print 'n initial queries: %s n kept: %s' % (len(ordered_queries), len(ordered_queries_filtered))
    for lf in t1.iter_leaves():
        #lf.add_face(AttrFace("name", fsize=20), 0, position="branch-right")
        lf.branch_vertical_margin = 0
        #data = [random.randint(0,2) for x in xrange(3)]

        for col, value in enumerate(ordered_queries_filtered):

            if head:
                if show_identity_values:
                    #'first row, print gene names'
                    #print 'ok!'
                    n = TextFace(' %s ' % str(value))
                    n.margin_top = 2
                    n.margin_right = 2
                    n.margin_left = 2
                    n.margin_bottom = 2
                    n.rotation = 270
                    n.vt_align = 2
                    n.hz_align = 2
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    #lf.add_face(n, col, position="aligned")
                    tss.aligned_header.add_face(n, col)
                else:
                    n = TextFace(' %s ' % str(value), fsize=6)
                    n.margin_top = 0
                    n.margin_right = 0
                    n.margin_left = 0
                    n.margin_bottom = 0
                    n.rotation = 270
                    n.vt_align = 2
                    n.hz_align = 2
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    #lf.add_face(n, col, position="aligned")
                    tss.aligned_header.add_face(n, col)
            try:

                identity_value = blast2data[lf.name][value][0]
                #print 'identity', lf.name, value, identity_value

                if lf.name != reference_accession:
                    if not accession2hit_filter:
                        color = rgb2hex(m_red.to_rgba(float(identity_value)))
                    else:
                        # if filter, color hits that are not in the filter in green
                        if value in accession2hit_filter[lf.name]:
                            color = rgb2hex(m_red.to_rgba(float(identity_value)))

                        else:
                            color = rgb2hex(m_green.to_rgba(float(identity_value)))
                else:
                    # reference taxon, blue scale
                    color = rgb2hex(m_blue.to_rgba(float(identity_value)))
                if not show_identity_values:
                    color = rgb2hex(m_blue.to_rgba(float(identity_value)))

            except:
                identity_value = 0
                color = "white"
            #print id_cutoff, float(identity_value)
            if show_identity_values:
                if float(identity_value) > float(id_cutoff):


                    if str(identity_value) == '100.00' or str(identity_value) == '100.0':
                        identity_value = '100'
                        n = TextFace("%s   " % identity_value)
                    else:
                        #    identity_value = str(round(float(identity_value), 1))

                        n = TextFace("%.2f" % round(float(identity_value), 2))
                    if float(identity_value) >95:
                        n.fgcolor = "white"

                    n.opacity = 1.
                else:
                    #print 'identity not ok', lf.name, identity_value, value
                    #print blast2data[lf.name]
                    identity_value = '-'
                    n = TextFace(' %s ' % str(identity_value))
                    n.opacity = 1.
                n.margin_top = 2
                n.margin_right = 2
                n.margin_left = 2
                n.margin_bottom = 2
                n.inner_background.color = color
                lf.add_face(n, col, position="aligned")
            else:
                # don't show identity values
                n = TextFace('  ')

                n.margin_top = 0
                n.margin_right = 0
                n.margin_left = 0
                n.margin_bottom = 0
                #n.color = color
                n.inner_background.color = color
                lf.add_face(n, col, position="aligned")


        try:
            lf.name = ' %s (%s)' % (id2description[lf.name], id2mlst[lf.name])
        except KeyError:
            lf.name = ' %s (%s)' % (lf.name, id2mlst[lf.name])
        head = False


    for n in t1.traverse():
       nstyle = NodeStyle()
       if n.support < 0.9:
           #mundo = TextFace("%s" % str(n.support))
           #n.add_face(mundo, column=1, position="branch-bottom")
           nstyle["fgcolor"] = "blue"
           nstyle["size"] = 6
           n.set_style(nstyle)
       else:
           nstyle["fgcolor"] = "red"
           nstyle["size"] = 0
           n.set_style(nstyle)


    print 'rendering tree'
    t1.render("profile.svg", dpi=1000, h=400, tree_style=tss)

    print blast2data
    print blast_result_file_list


def get_accession_filter_from_blast_list(tab_blast_results, identity_cutoff=80):
    accession2locus = {}
    for one_blast in tab_blast_results:
        with open(one_blast, 'r') as f:
            for row in f:
                data = row.rstrip().split()
                if float(data[2])>=identity_cutoff:
                    if data[1] not in accession2locus:
                        accession2locus[data[1]] = [data[0]]
                    else:
                        accession2locus[data[1]].append(data[0])

    return accession2locus



def main(input_reference,
         input_queries_folder,
         blast_file,
         mlst_scheme,
         input_gbk,
         skip_parsnp=False,
         skip_blast=False,
         input_tree=None,
         skip_mlst=False,
         check_overlap=False,
         id_cutoff=80,
         blast_type='tblastn',
         reference_accession='-',
         blast_filter=False,
         show_identity_values=True,
         accession2description=False):

    import shell_command
    import os
    import sys
    import glob
    from Bio import SeqIO
    sys.stdout.write('Building tree using parsnp...\n')

    wd = os.getcwd()

    ordered_queries = [record.name for record in SeqIO.parse(blast_file, 'fasta')]

    #print 'input fasta folder', input_queries_folder
    fasta_folder = os.path.abspath(input_queries_folder)
    reference_file = os.path.abspath(input_reference)

    if blast_filter:
        accession2hit_filter=get_accession_filter_from_blast_list(blast_filter)
    else:
        accession2hit_filter=False

    #print 'fasta folder', fasta_folder
    pp = fasta_folder + '/*fna'
    #print pp
    fasta_files = glob.glob(pp)
    if len(fasta_files) == 0:
        pp = fasta_folder + '/*ffn'
        #print pp
        fasta_files = glob.glob(pp)
        if len(fasta_files) == 0:
            raise ('could not find fasta files')
    #print 'fasta files', fasta_files

    reference_phylogeny_folder = os.path.join(wd, 'reference_parsnp_phylogeny')
    reference_phylogeny = os.path.join(reference_phylogeny_folder, 'parsnp_edit.tree')
    out_mlst = os.path.join(wd, 'mlst_results/mlst.tab')
    if not skip_parsnp and not input_tree:
        cmd = 'parsnp -r %s -d %s -p 6 -c -o parsnp_tree' % (reference_file, fasta_folder)
        #print cmd
        out, err, code = shell_command.shell_command(cmd)
        print code, err
        print out

        if not os.path.exists(reference_phylogeny_folder):
            os.mkdir(reference_phylogeny_folder)
        sys.stdout.write('Editing parsnp phylogeny...\n')

        parsnp_raw_phylogeny = os.path.join(wd, 'parsnp_tree/parsnp.tree')
        cmd = 'cat %s | sed  "s/.fna//g"> %s' % (parsnp_raw_phylogeny, reference_phylogeny)
        out, err, code = shell_command.shell_command(cmd)
        #print code, err
        cmd = '''sed -i "s/\'//g" %s''' % reference_phylogeny
        out, err, code = shell_command.shell_command(cmd)
        #print cmd, '##################'
        #print code, err

        cmd = '''sed -i "s/.ref//" %s''' % (reference_phylogeny)
        out, err, code = shell_command.shell_command(cmd)
        #print code, err

    os.chdir(wd)
    if not skip_mlst:
        sys.stdout.write('Identifying mlst...\n')
        if not os.path.exists(os.path.join(wd, 'mlst_results')):
            os.mkdir(os.path.join(wd, 'mlst_results'))
        all_fasta = ' '.join(fasta_files)

        cmd = 'mlst --quiet --nopath --scheme %s %s > %s' % (mlst_scheme, all_fasta, out_mlst)
        out, err, code = shell_command.shell_command(cmd)
        cmd = 'sed -i "s/.fna//g" %s' % out_mlst
        out, err, code = shell_command.shell_command(cmd)
        cmd = 'sed -i "s/.ffn//g" %s' % out_mlst
        out, err, code = shell_command.shell_command(cmd)

    #print err, code
    if input_tree:
        reference_phylogeny = input_tree



    sys.stdout.write('Blasting %s...\n' % blast_file)

    blast_file_fullpath = os.path.join(wd, os.path.basename(blast_file))

    blastp_folder = os.path.join(wd, 'blastp_results')
    if not os.path.exists(blastp_folder):
        os.mkdir(blastp_folder)
    os.chdir(input_queries_folder)

    blast_best_hit_results = []
    for genome_file in fasta_files:
        out_file = 'blast_'+ os.path.basename(genome_file).split('.')[0]
        outpath = os.path.join(blastp_folder, out_file + '.tab')
        best_hit_path = os.path.join(blastp_folder, 'uniq_' + out_file + '.tab')
        blast_best_hit_results.append(best_hit_path)
        if not skip_blast:
            cmd = 'formatdb -i %s -p F' % genome_file
            out, err, code = shell_command.shell_command(cmd)
            if blast_type == 'tblastn':
                from Bio.Blast.Applications import NcbitblastnCommandline

                #print err, code
                out, err, code = shell_command.shell_command('export BLASTDB=$BLASTDB:%s' % blastp_folder)


                blastp_cline = NcbitblastnCommandline(query= blast_file_fullpath,
                                                    db=genome_file,
                                                    evalue=0.001,
                                                    outfmt=6,
                                                    out=outpath)

            elif blast_type == 'blastn':
                from Bio.Blast.Applications import NcbiblastnCommandline
                blastp_cline = NcbiblastnCommandline(query= blast_file_fullpath,
                                                    db=genome_file,
                                                    evalue=0.001,
                                                    outfmt=6,
                                                    out=outpath)

            else:
                raise ('unsupported blast type')

            stdout, stderr = blastp_cline()
            result_handle = open(outpath, 'r')
            best_hit_handle = open(best_hit_path, 'w')

            hit_list = []
            for line in result_handle:
                if line.split('\t')[0] in hit_list:
                    continue
                else:
                    hit_list.append(line.split('\t')[0])
                    best_hit_handle.write(line)
            best_hit_handle.close()
    os.chdir(wd)


    if accession2description:
        id2description = {}
        with open(accession2description, 'r') as f:
            for row in f:
                data = row.rstrip().split('\t')
                id2description[data[0]] = data[1]
    else:
        id2description = gbk2accessiontodefinition.get_coressp(input_gbk)
    #else:
    #    IOError('either provide id2description og gbk files')

    accession2st_type = parse_mlst_results(out_mlst)

    #print reference_phylogeny

    plot_blast_result(reference_phylogeny,
                      blast_best_hit_results,
                      id2description,
                      accession2st_type,
                      check_overlap,
                      ordered_queries,
                      id_cutoff,
                      reference_accession=reference_accession,
                      accession2hit_filter=accession2hit_filter,
                      show_identity_values=show_identity_values)


if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb
    import gbk2accessiontodefinition
    parser = argparse.ArgumentParser()
    parser.add_argument("-r",'--reference', type=str, help="reference genome")
    parser.add_argument("-f",'--input_queries_folder', type=str,help="genomes folder")
    parser.add_argument("-b",'--blast_file', type=str, help="fasta file with sequences to blast and map on the phylogeny")
    parser.add_argument("-g",'--input_gbk', type=str, help="corresponding gbk files to get id correspondance", nargs="+")
    parser.add_argument("-m",'--mlst_shemes', type=str, help="mlst sheme /t seeman script/")
    parser.add_argument("-s",'--skip_parsnp', action='store_true', help="skip parsnp part")
    parser.add_argument("-l",'--skip_blast', action='store_true', help="skip blast part")
    parser.add_argument("-t",'--input_tree', type=str, help="user provided tree /no parsnp tree/", default=None)
    parser.add_argument("-n",'--skip_mlst', action='store_true', help="skip mlst part")
    parser.add_argument("-c",'--id_cutoff', default=80, help="identity cutoff /default 80 percent/", type=int)
    parser.add_argument("-bt",'--blast_type', default='tblastn', help="blast type (default: tblastn. Alternative: blastn)")
    parser.add_argument("-re",'--reference_accession', default='-', help="reference accession (colored in blue scale)")
    parser.add_argument("-bf",'--blast_filter', default=False, help="blast filter (color in green hits found in the provided blast files)", nargs='+')
    parser.add_argument("-id",'--dont_display_id', action="store_false", help="don't show identity values", default=True)
    parser.add_argument("-a",'--accession2description',default=False, type=str,help="tab file with accessions and corresponding descriptions for leaf labels")


    args = parser.parse_args()

    if 'ref' in args.reference:
        import shutil
        shutil.copyfile(args.reference, 'temp.fna')
        ref = 'temp.fna'
    else:
        ref = args.reference




    main(ref,
         args.input_queries_folder,
         args.blast_file,
         args.mlst_shemes,
         args.input_gbk,
         args.skip_parsnp,
         args.skip_blast,
         args.input_tree,
         args.skip_mlst,
         False,
         float(args.id_cutoff),
         blast_type=args.blast_type,
         reference_accession=args.reference_accession,
         blast_filter=args.blast_filter,
         show_identity_values=args.dont_display_id,
         accession2description=args.accession2description)



    #id2description = gbk2accessiontodefinition.get_coressp(args.input_gbk)
    #accession2st_type = parse_mlst_results(args.st_file)
    #plot_blast_result(args.tree, args.blast_files, id2description, accession2st_type)
