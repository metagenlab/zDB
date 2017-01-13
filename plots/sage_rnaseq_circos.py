#!/usr/bin/env python
# python 2.7.5 requires biopython


def circos_rnaseq_abundance():
    a = open('/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/count_filter.tab', 'r')
    b = open('/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/locu_start_stop.tab', 'r')

    locus2start_stop = {}
    locus2count = {}

    for i, line in enumerate(a):
        if i == 0:
            data = line.rstrip().split('\t')[1:len(line)]
            samples = [i[1:-1] for i in data]
        else:
            data = line.rstrip().split('\t')
            locus2count[data[0][1:-1]] = {}
            for y, sample in enumerate(samples):
                locus2count[data[0][1:-1]][sample] = data[y+1]
    for line in b:
        data = line.rstrip().split('\t')
        locus2start_stop[data[-1]] = data[0:-1]
    #print locus2count
    print locus2start_stop

    cond2file = {}
    for sample in samples:
        cond2file[sample] = open('/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/circos_%s.tab' % sample, 'w')

    contig = "S5_genome"
    for locus in locus2count:
        for sample in samples:
            cond2file[sample].write("%s\t%s\t%s\t%s\tid=%s\n" % (contig,
                                                        locus2start_stop[locus][0],
                                                        locus2start_stop[locus][1],
                                                        locus2count[locus][sample],
                                                                 locus2start_stop[locus][3]))



def circos_locus2heatmap(expression_data, locus2start_stop, outfile):

    a = open('/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/count_filter.tab', 'r')
    b = open(locus2start_stop, 'r')
    c = open(expression_data, 'r')

    locus2start_stop = {}
    locus2count = {}
    locus2fold_change={}

    for i, line in enumerate(a):
        if i == 0:
            data = line.rstrip().split('\t')[1:len(line)]
            samples = [i[1:-1] for i in data]
        else:
            data = line.rstrip().split('\t')
            locus2count[data[0][1:-1]] = {}
            for y, sample in enumerate(samples):
                locus2count[data[0][1:-1]][sample] = data[y+1]
    for line in b:
        data = line.rstrip().split('\t')
        locus2start_stop[data[-1]] = data[0:-1]

    locus_list = []
    for line in c:
        data = line.rstrip().split('\t')
        locus = data[1][1:-1]
        locus_list.append(locus)
        locus2fold_change[locus] = data[2]


    #print locus2count
    #print locus2start_stop

    f = open(outfile, 'w')




    contig = "S5_genome"
    for locus in locus2start_stop:
        if locus2start_stop[locus][0] == "start":
            continue
        if locus in locus_list:

            f.write("%s\t%s\t%s\t%s\tid=%s\n" % (contig,
                                                locus2start_stop[locus][0],
                                                locus2start_stop[locus][1],
                                                locus2fold_change[locus],
                                                locus2start_stop[locus][3]))
        else:
            f.write("%s\t%s\t%s\t%s\tid=%s\n" % (contig,
                                                locus2start_stop[locus][0],
                                                locus2start_stop[locus][1],
                                                "0",
                                                locus2start_stop[locus][3]))


def circos_locus2hilight(locus_data, locus2start_stop, outfile):
    a = open('/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/count_filter.tab', 'r')
    b = open(locus2start_stop, 'r')

    locus2start_stop = {}
    locus2count = {}

    for i, line in enumerate(a):
        if i == 0:
            data = line.rstrip().split('\t')[1:len(line)]
            samples = [i[1:-1] for i in data]
        else:
            data = line.rstrip().split('\t')
            locus2count[data[0][1:-1]] = {}
            for y, sample in enumerate(samples):
                locus2count[data[0][1:-1]][sample] = data[y+1]
    for line in b:
        data = line.rstrip().split('\t')
        locus2start_stop[data[-1]] = data[0:-1]
    #print locus2count
    #print locus2start_stop

    f = open(outfile, 'w')

    contig = "S5_genome"
    for locus in locus_data:
        print locus
        f.write("%s\t%s\t%s\tid=%s\n" % (contig,
                                         locus2start_stop[locus][0],
                                         locus2start_stop[locus][1],
                                                 locus2start_stop[locus][3]))


#circos_rnaseq_abundance()


'''
locus_list = [line.rstrip().split('\t')[1] for line in open("/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/LM_WL_top200.txt", 'r')]
locus_list2 = [line.rstrip().split('\t')[1] for line in open("/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/LM_WR_top200.txt", 'r')]
locus_list3 = [line.rstrip().split('\t')[1] for line in open("/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/LM_SA_top200.txt", 'r')]
locus_list4 = [line.rstrip().split('\t')[1] for line in open("/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/WL_WR_top200.txt", 'r')]



circos_locus2highlight(locus_list, '/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/locu_start_stop.tab',
                       '/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/circos_LM_WL_top_200.tab')
circos_locus2highlight(locus_list2, '/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/locu_start_stop.tab',
                       '/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/circos_LM_WR_top_200.tab')

circos_locus2highlight(locus_list3, '/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/locu_start_stop.tab',
                       '/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/circos_LM_SA_top_200.tab')

circos_locus2highlight(locus_list4, '/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/locu_start_stop.tab',
                       '/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/circos_WL_WR_top_200.tab')
'''



circos_locus2heatmap("/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/wl_rw_all.tab",
                     '/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/locu_start_stop.tab',
                     '/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/circos_WL_WR_fold_change.tab')


circos_locus2heatmap("/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/lm_wr_all.tab",
                     '/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/locu_start_stop.tab',
                     '/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/circos_LM_WR_fold_change.tab')


circos_locus2heatmap("/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/lm_wl_all.tab",
                     '/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/locu_start_stop.tab',
                     '/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/circos_LM_WL_fold_change.tab')


circos_locus2heatmap("/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/lm_sa_all.tab",
                     '/home/trestan/ownCloud/SAGE/2015/rna_seq/genome/locu_start_stop.tab',
                     '/home/trestan/ownCloud/SAGE/2015/rna_seq/comparative_analysis/gbk2circos/circos_LM_SA_fold_change.tab')

