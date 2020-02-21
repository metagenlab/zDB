#!/usr/bin/python

import rpy2.robjects.numpy2ri

import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr
import rpy2.rlike.container as rlc
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import pandas.rpy.common as com


def accession2dna_seg_genome(accession, biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    import numpy
    import pandas

    server, db = manipulate_biosqldb.load_db(biodb)
    sql = 'select locus_tag, start, stop, strand from orthology_detail where accession="%s"' % (biodb, accession)

    data = numpy.array([list(i) for i in server.adaptor.execute_and_fetchall(sql,)])
    columns = ['name', 'start', 'end', 'strand']

    df = pandas.DataFrame(data, columns=columns)
    return df

def get_genome_seg(accession, biodb):

    from chlamdb.biosqldb import manipulate_biosqldb
    import pandas

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select genome_size from genomes_info where accession="%s";' % (biodb, accession)

    genome_size = server.adaptor.execute_and_fetchall(sql,)[0][0]
    data = [(accession, 0, genome_size, '-')]
    columns = ['name', 'start', 'end', 'strand']

    df = pandas.DataFrame(data, columns=columns)
    return df


def get_pairwise_connexions(accession_1, accession_2, biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    import numpy
    import pandas

    server, db = manipulate_biosqldb.load_db(biodb)

    sql1 = 'select seqfeature_id, start, stop from orthology_detail where accession in ("%s","%s") ' % (biodb,
                                                                                                                   accession_1,
                                                                                                                   accession_2)
    seqfeature_id2location = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql1,))

    print seqfeature_id2location.keys()[0:10]

    sql2 = 'select accession, taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
           ' where t1.name="%s"' % biodb
    print sql2


    accession2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    comp1_sql = 'select locus_1,locus_2,identity from (select * from ' \
                ' comparative_tables_identity_closest_homolog2 where taxon_1=%s and taxon_2=%s) A ' \
                ' inner join orthology_detail B on A.locus_1=B.seqfeature_id;' % (biodb,
                                                                                              accession2taxon_id[accession_1],
                                                                                              accession2taxon_id[accession_2],
                                                                                              biodb)


    data = server.adaptor.execute_and_fetchall(comp1_sql,)

    comparison_table = []
    for row in data:
        print row
        try:
            start1 = seqfeature_id2location[int(row[0])][0]
            stop1 = seqfeature_id2location[int(row[0])][0]
            start2 = seqfeature_id2location[int(row[1])][0]
            stop2 = seqfeature_id2location[int(row[1])][0]
            identity = row[2]
            comparison_table.append([start1, stop1, start2, stop2, identity])
        except:
            pass
    data = numpy.array(comparison_table)
    columns = ['start1', 'end1', 'start2', 'end2', 'identity']
    df = pandas.DataFrame(data, columns=columns)
    return df

def accession_list2plot(accession_list, biodb):

    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    import rpy2.robjects.numpy2ri as numpy2ri


    robjects.r('''
        dna_seg_dataframes <- list()
        dna_comparison_dataframes <- list()
    ''')

    print 'getting table 1'
    seg_table_list = []
    for i, accession in enumerate(accession_list):
        #table1 = taxon2dna_seg_genome(accession_list[0], biodb)
        table1 = accession2dna_seg_genome(accession, biodb)#get_genome_seg(accession, biodb)
        seg_table_list.append(com.convert_to_r_dataframe(table1))

    robjects.r.assign('dna_seg_dataframes', seg_table_list)
    comparison_tables = []
    for x in range(0, len(accession_list)-1):

        first_accession = accession_list[x]
        second_accession = accession_list[x+1]

        #table2 = taxon2dna_seg_genome(accession_list[1], biodb)
        #table2 = get_genome_seg(accession_list[1], biodb)

        print 'getting table 3'
        table_comparison = get_pairwise_connexions(first_accession, second_accession, biodb)
        comparison_tables.append(com.convert_to_r_dataframe(table_comparison))


    robjects.r.assign('dna_comparison_dataframes', comparison_tables)

    #robjects.r.assign('df3', com.convert_to_r_dataframe(table_comparison))

    robjects.r('''
        library(genoPlotR)
        library(Cairo)

        print (length(dna_seg_dataframes))
        print (dna_seg_dataframes)
        print (length(dna_comparison_dataframes))
        print (dna_comparison_dataframes)

        for (i in 1:length(dna_seg_dataframes)){
        print(i)
            dna_seg_dataframes[[i]]$start <- as.numeric(as.character(dna_seg_dataframes[[i]]$start))
            dna_seg_dataframes[[i]]$end <- as.numeric(as.character(dna_seg_dataframes[[i]]$end))
            dna_seg_dataframes[[i]] <- dna_seg(dna_seg_dataframes[[i]])
            dna_seg_dataframes[[i]]$gene_type <- rep("side_blocks", length(dna_seg_dataframes[[i]]$start))
        }

        for (i in 1:length(dna_comparison_dataframes)){
            print(i)
            dna_comparison_dataframes[[i]]$start1 <- as.numeric(as.character(dna_comparison_dataframes[[i]]$start1))
            dna_comparison_dataframes[[i]]$end1 <- as.numeric(as.character(dna_comparison_dataframes[[i]]$end1))
            dna_comparison_dataframes[[i]]$start2 <- as.numeric(as.character(dna_comparison_dataframes[[i]]$start2))
            dna_comparison_dataframes[[i]]$end2 <- as.numeric(as.character(dna_comparison_dataframes[[i]]$end2))
            dna_comparison_dataframes[[i]] <- as.comparison(dna_comparison_dataframes[[i]])

        }
        #df1$start <- as.numeric(as.character(df1$start))
        #df1$end <- as.numeric(as.character(df1$end))
        #df2$start <- as.numeric(as.character(df2$start))
        #df2$end <- as.numeric(as.character(df2$end))
        #df3$start1 <- as.numeric(as.character(df3$start1))
        #df3$end1 <- as.numeric(as.character(df3$end1))
        #df3$start2 <- as.numeric(as.character(df3$start2))
        #df3$end2 <- as.numeric(as.character(df3$end2))



        #print(head(df1))
        #print(head(df2))
        #print(head(df3))

        #dna_seg1 <- dna_seg(df1)
        #dna_seg2 <- dna_seg(df2)
        #comparison <- as.comparison(df3)
        #dna_segs <- list(dna_seg1,dna_seg2)

        CairoPDF('test2.pdf',height=6,width=15)
        plot_gene_map(dna_segs=dna_seg_dataframes, comparisons=dna_comparison_dataframes)
        dev.off()
    ''')





if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--seq_proj_GOLD', default=False, type=str, help="GOLD sequencing project ID")

    args = parser.parse_args()

    accession_list2plot(['NZ_LN879502','NC_005861', 'NC_015702', 'NC_014225', 'NC_015713'],'2017_06_29b_motile_chlamydiae')
