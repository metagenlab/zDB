#!/usr/bin/env python

import rpy2.robjects.numpy2ri
import rpy2.robjects as robjects
rpy2.robjects.numpy2ri.activate()





def gc_coverage_plot(samtool_depth_file, contigs_file, blast_file=False, column1=1, column2=2, main=False, highlight=False):

    import os
    if not main:
        main = os.path.basename(samtool_depth_file)

    import shell_command

    out, err, code = shell_command.shell_command("infoseq -auto -only -Name -length -pgc %s > /tmp/gc.tab" % contigs_file)

    print out
    print err
    print code


    if highlight:
        highlight_code = """
         gc_coverage_table$color <- rep(rgb(1, 0, 0,0.5), length(gc_coverage_table[,1]))
        highlight_table <- read.table("%s", header=FALSE)
        m <- match(highlight_table[,1], gc_coverage_table$Name)
        gc_coverage_subset <- gc_coverage_table[m,]
        print("subset")
        print(m)
        gc_coverage_table[m,]$color<-rgb(0, 0, 1,0.5)

        """ % highlight

        highlight_code2 = """

        m <- match(highlight_table[,1], gc_coverage_table_2m$Name)
        print("subset m2")
        print(m)
        gc_coverage_subset2 <- gc_coverage_table_2m[m,]

        """

    else:
        highlight_code = ''
        highlight_code2 = ''

    print 'high', highlight

    if not blast_file:
        robjects.r("""

        #library(Cairo)
        library(R.utils)

        if (isGzipped("%s")){
            print('Gzipped file')
            all_depth <- read.table(gzfile('%s'), header=FALSE)
        }else{
            print('Not Gzipped')
            all_depth <- read.table('%s', header=FALSE)
        }


        contigs_depth<- aggregate(all_depth["V3"],b=all_depth[c("V1")],FUN=median)
        contigs_gc <- read.table("/tmp/gc.tab", header=TRUE)

        gc_coverage_table <-cbind(contigs_gc,coverage=contigs_depth[match(contigs_gc$Name, contigs_depth$V1),2])
        w<-which(gc_coverage_table$Length >=1000)
        gc_coverage_table <- gc_coverage_table[w,]

        

        write.table(gc_coverage_table, 'gc_coverage_table.tab', sep="\t", row.names=F)


    %s

     svg("gc_cov_buble.svg", width = 12, height = 12,)
    symbols(x=gc_coverage_table[,3], y= gc_coverage_table[,4], circles=gc_coverage_table[,2], inches=1/3, ann=F, bg=rgb(1, 0, 0,0.5), fg=rgb(1, 0, 0,0.5), main=%s)
    if (any("gc_coverage_subset" %%in%% ls())) {

     symbols(x=gc_coverage_table[,3], y= gc_coverage_table[,4], circles=gc_coverage_table[,2], inches=1/3, ann=F, bg=gc_coverage_table$color, fg=gc_coverage_table$color, add = TRUE)
l <- gsub('(^[^_]+_[^_]+)_(.*)$', '\\\\1', gc_coverage_subset$Name)
     text(x=gc_coverage_subset[,3], y=gc_coverage_subset[,4], labels = l)
     }else{

     print ('a') }

     dev.off()

     cov_biggest <- gc_coverage_table[which(gc_coverage_table$Length==max(gc_coverage_table$Length)),4]
     print('cov biggest:')
     print(cov_biggest)
     w <- which(gc_coverage_table[,4]< (4*cov_biggest))
     gc_coverage_table_2m <- gc_coverage_table[w,]


     %s

     svg("gc_cov_buble_2m.svg", width = 12, height = 12,)
     symbols(x=gc_coverage_table_2m[,3], y= gc_coverage_table_2m[,4], circles=gc_coverage_table_2m[,2], inches=1/3, ann=F, bg=rgb(1, 0, 0,0.5), fg=rgb(1, 0, 0,0.5), main=%s)



    if (any("gc_coverage_subset" %%in%% ls())) {



     symbols(x=gc_coverage_table_2m[,3], y= gc_coverage_table_2m[,4], circles=gc_coverage_table_2m[,2], inches=1/3, ann=F, bg=gc_coverage_table_2m$color, fg=gc_coverage_table_2m$color, add = TRUE)
     l <- gsub('(^[^_]+_[^_]+)_(.*)$', '\\\\1', gc_coverage_subset2$Name)
     print (l)
     text(x=gc_coverage_subset2[,3], y=gc_coverage_subset2[,4], labels = l)
     }else{

     print ('a') }

     dev.off()



                   """ % (samtool_depth_file,
                          samtool_depth_file,
                          samtool_depth_file,
                          highlight_code,
                          main,
                          highlight_code2,
                          main))
    else:

        robjects.r("""

        #library(Cairo)
        library(R.utils)

        if (isGzipped("%s")){
            print('Gzipped file')
            all_depth <- read.table(gzfile('%s'), header=FALSE)
        }else{
            print('Not Gzipped')
            all_depth <- read.table('%s', header=FALSE)
        }

        blast_file <- read.table("%s", header=FALSE, sep="\t")[,c(2,6)]
        contigs_depth<- aggregate(all_depth["V3"],b=all_depth[c("V1")],FUN=median)
        contigs_gc <- read.table("/tmp/gc.tab", header=TRUE)

        gc_coverage_table <-cbind(contigs_gc,coverage=contigs_depth[match(contigs_gc$Name, contigs_depth$V1),2])
        w<-which(gc_coverage_table$Length >=1000)
        gc_coverage_table <- gc_coverage_table[w,]

        gc_coverage_table$taxon <- blast_file[,2][match(gc_coverage_table$Name, blast_file[,1])]
        print (is.na(gc_coverage_table$taxon))
        gc_coverage_table$taxon <- as.character(gc_coverage_table$taxon)
        gc_coverage_table$taxon[is.na(gc_coverage_table$taxon)] <- 'undefined'
        gc_coverage_table$taxon <- as.factor(gc_coverage_table$taxon)
        
        write.table(gc_coverage_table, 'gc_coverage_table.tab', sep="\t", row.names=F)

     svg("gc_cov_buble.svg", width = 12, height = 12,)
     symbols(x=gc_coverage_table[,3], y= gc_coverage_table[,4], circles=gc_coverage_table[,2], inches=1/3, ann=F, bg=gc_coverage_table$taxon, fg=gc_coverage_table$taxon, main=%s)
     dev.off()

     cov_biggest <- gc_coverage_table[which(gc_coverage_table$Length==max(gc_coverage_table$Length)),4]
     print('cov biggest:')
     print(cov_biggest)
     w <- which(gc_coverage_table[,4]< (4*cov_biggest))
     gc_coverage_table_2m <- gc_coverage_table[w,]

     svg("gc_cov_buble_2m.svg", width = 12, height = 12,)
     symbols(x=gc_coverage_table_2m[,3], y= gc_coverage_table_2m[,4], circles=gc_coverage_table_2m[,2], inches=1/3, ann=F, bg=gc_coverage_table$taxon, fg=gc_coverage_table$taxon, main=%s)
     dev.off()



                   """ % (samtool_depth_file,
                          samtool_depth_file,
                          samtool_depth_file,
                          blast_file,
                          main,
                          main))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_samtools_depth', type=str, help="input samtools depth files")
    parser.add_argument("-m", '--input_contigs', type=str, help="input contigs")
    parser.add_argument("-b", '--blast_file', type=str, help="file with blast columns")
    parser.add_argument("-1", '--column1', type=str, help="contig names column", default=2)
    parser.add_argument("-2", '--column2', type=str, help="classification column", default=6)
    parser.add_argument("-l", '--highlight', type=str, help="highlight some contigs listed in input file", default=False)

    args = parser.parse_args()

    gc_coverage_plot(args.input_samtools_depth, args.input_contigs, args.blast_file, args.column1, args.column2, False, args.highlight)
