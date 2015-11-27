#!/usr/bin/env python

import rpy2.robjects.numpy2ri
import rpy2.robjects as robjects
rpy2.robjects.numpy2ri.activate()





def gc_coverage_plot(samtool_depth_file, contigs_file, main=False):

    if not main:
        main = samtool_depth_file

    import shell_command

    out, err, code = shell_command.shell_command("infoseq -auto -only -Name -length -pgc %s > /tmp/gc.tab" % contigs_file)

    print out


    robjects.r("""

    library(Cairo)
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
    w<-which(gc_coverage_table$Length >1000)
    gc_coverage_table <- gc_coverage_table[w,]

    print(head(gc_coverage_table))

 CairoSVG("gc_cov_buble.svg", width = 12, height = 12,)
symbols(x=gc_coverage_table[,3], y= gc_coverage_table[,4], circles=gc_coverage_table[,2], inches=1/3, ann=F, bg="steelblue2", fg="red", main=%s)
 dev.off()

 cov_biggest <- gc_coverage_table[which(gc_coverage_table$Length==max(gc_coverage_table$Length)),4]
 print('cov biggest:')
 print(cov_biggest)
 w <- which(gc_coverage_table[,4]< (4*cov_biggest))
 gc_coverage_table_2m <- gc_coverage_table[w,]

 CairoSVG("gc_cov_buble_2m.svg", width = 12, height = 12,)
symbols(x=gc_coverage_table_2m[,3], y= gc_coverage_table_2m[,4], circles=gc_coverage_table_2m[,2], inches=1/3, ann=F, bg="steelblue2", fg="red", main=%s)
 dev.off()



               """ % (samtool_depth_file,
                      samtool_depth_file,
                      samtool_depth_file,
                      main,
                      main))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_samtools_depth', type=str, help="input samtools depth files")
    parser.add_argument("-m", '--input_contigs', type=str, help="input contigs")


    args = parser.parse_args()

    gc_coverage_plot(args.input_samtools_depth, args.input_contigs)
