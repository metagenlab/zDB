#!/usr/bin/env python

import numpy
import rpy2.robjects.numpy2ri
import rpy2.robjects as robjects
rpy2.robjects.numpy2ri.activate()





def samtools_depth2coverage_plot(samtool_depth_file, main=False):

    if not main:
        import os
        main = os.path.basename(samtool_depth_file)

    robjects.r("""

    #library(Cairo)
    library(R.utils)


    if (isGzipped("%s")){
        print('Gzipped file')
        my_file <- read.table(gzfile('%s'), header=FALSE)
    }else{
        my_file <- read.table('%s', header=FALSE)
    }

    contig_lengths <- table(factor(my_file$V1, levels=unique(my_file$V1)))

    contig_limits = c()
    length = 0
    for (i in contig_lengths){
        length <- length + as.numeric(i)
        print('length')
        print(length)
        contig_limits <- c(contig_limits,length)
    }
    print('ok')
    print(contig_limits)
    cov_data <- my_file$V3

    median_depth = median(cov_data)
    max_depth = max(cov_data)

    newlength <- ceiling(length(cov_data)/100)*100
    cov_data[newlength] <- NA
    cov_matrix <- matrix(cov_data,nrow=100)
    cov_100bp <- colMeans(cov_matrix, na.rm=TRUE)
    write.table(as.data.frame(cov_100bp), "coverage_100bp.tab", sep="\t")
    pdf('%s_coverage.pdf')
        #par(mfrow=c(1,2))
        x<-seq(100,newlength,100)
        print(max_depth)
        plot(x, cov_100bp,type='l',col='light grey', las=2, main='%s', ylim=c(0,max_depth), xlab="Position", ylab="Sequencing depth")
        f1001 <- rep(1/1001,1001)
        y_sym <- filter(cov_100bp, f1001,sides=2)
        lines(x,y_sym,col="blue")
        text(length(cov_data)*0.2,max_depth*0.2, paste("Median depth:", median_depth), col="blue")
        abline(v=contig_limits, col=rgb(1, 0, 0, 0.5), lty=3, lwd=0.5)
        abline(h=10, col="red")
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
    parser.add_argument("-m", '--main', type=str, help="plot title", default=False)


    args = parser.parse_args()

    samtools_depth2coverage_plot(args.input_samtools_depth, args.main)
