#!/usr/bin/env python

import numpy
import rpy2.robjects.numpy2ri
import rpy2.robjects as robjects


def samtools_depth2coverage_plot(samtool_depth_file, fasta_file=False, main=False):

    if not main:
        import os
        main = os.path.basename(samtool_depth_file)

    if fasta_file:
        print ('fasta file --------------------------------')
        from Bio import SeqIO
        concatenated_seq = ''
        with open(fasta_file, 'r') as f:
            records = SeqIO.parse(f, 'fasta')
            for record in records:
                concatenated_seq+=str(record.seq)
        robjects.r.assign('dna_sequence', concatenated_seq)
        robjects.r.assign('show_GC', True)
    else:
        robjects.r.assign('show_GC', False)

    robjects.r("""

    #library(Cairo)
    library(R.utils)
    library(grid)
    library(gridBase)
    library(seqinr)
    library(pracma)
    library(zoo)
    
    slidingwindowplot_GC <- function(windowsize, windowsize2, inputseq, contig_limits)
    {
       starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
       n <- length(starts)    # Find the length of the vector "starts"
       chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
       for (i in 1:n) {
            chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
            chunkGC <- GC(chunk)*100
            chunkGCs[i] <- chunkGC
       }
       plot(starts,chunkGCs,type="l", ylab="GC(%%)", xaxt='n', xaxs="i", frame.plot = FALSE, col="blue", lwd = 0.5)
       axis(side = 2)

       #roll <- runmed(chunkGCs, k=windowsize2)
       #lines(starts, roll, col="blue")
       
       abline(v=contig_limits, col=rgb(1, 0, 0, 0.5), lty=3, lwd=0.5)
       #axis(side=2, at=c(20,30,40,50,60), labels=c("20%%","30%%","40%%","50%%","60%%"))
       abline(h=GC(inputseq))
    }

    
    slidingwindowplot_depth <- function(windowsize1, windowsize2, depth, contig_limits)
    {
       starts1 <- seq(1, length(depth)-windowsize1, by = windowsize1)
       n1 <- length(starts1)    # Find the length of the vector "starts"
       chunks_mean_depth1 <- numeric(n1) # Make a vector of the same length as vector "starts", but just containing zeroes
       for (i in 1:n1) {
            chunk1 <- depth[starts1[i]:(starts1[i]+windowsize1-1)]
            chunk_mean_depth1 <- mean(chunk1)
            chunks_mean_depth1[i] <- chunk_mean_depth1
       }

        median_depth = median(chunks_mean_depth1)
        max_depth = max(chunks_mean_depth1)
        
        w1 <- which(chunks_mean_depth1>(3*median_depth))
        chunks_mean_depth1[w1] <- 3*median_depth     

        plot(starts1, chunks_mean_depth1,type='h',col='light grey', xaxs="i", las=2, main='', ylim=c(0,3*median_depth), xlab="Position", ylab="Sequencing depth", frame.plot = FALSE)
        axis(side = 2)

        roll <- runmed(depth, k=windowsize2)

        lines(roll, col="blue")
            
        #text(length(depth)*0.2,max_depth*0.2, paste("Median depth:", median_depth), col="blue")
       
        abline(v=contig_limits, col=rgb(1, 0, 0, 0.5), lty=3, lwd=0.5)
     }

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
        contig_limits <- c(contig_limits,length)
    }

    cov_data <- my_file$V3

    median_depth = median(cov_data)
    max_depth = max(cov_data)

    pdf('%s_coverage2.pdf', height=9, width=25)
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(2, 1, heights=unit(c(0.5,1), rep("null", 2)))))

            pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
            par(new=TRUE, fig=gridFIG(), las=1, mar=c(5, 5, 0, 5)) 

            slidingwindowplot_depth(150,5001, cov_data, contig_limits)

            upViewport()         

            if (show_GC != FALSE) {
                print("SHOW ---------------------------------------------------------------------------")
                    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
                    par(new=TRUE, fig=gridFIG(), las=1, mar=c(0, 5, 3, 5)) 
                    myseq <- unlist(strsplit(dna_sequence, ""))
                    slidingwindowplot_GC(2000, 21, myseq, contig_limits)
                    # ,xlab="Nucleotide start position",ylab="GC content"
                    #axis(1, at = seq(0, length(myseq), 10000), labels = seq(0, length(myseq), 10000)/1000, las=2)
                    upViewport()
            }

    dev.off()
               """ % (samtool_depth_file,
                      samtool_depth_file,
                      samtool_depth_file,
                      main))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_samtools_depth', type=str, help="input samtools depth files")
    parser.add_argument("-f", '--fasta_contig_file', type=str, help="fasta_contig_file")
    parser.add_argument("-m", '--main', type=str, help="plot title", default=False)


    args = parser.parse_args()

    samtools_depth2coverage_plot(args.input_samtools_depth, args.fasta_contig_file, args.main)
