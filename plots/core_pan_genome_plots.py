#!/usr/bin/python

import rpy2.robjects.numpy2ri

import numpy as np
from rpy2.robjects import pandas2ri
pandas2ri.activate()

def core_pan_genome_plot(numpy_matrix,
                          header="",
                          xlab="",
                          ylab="",
                          reverse=False,
                          output_path="~/test.svg"):

        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri
        rpy2.robjects.numpy2ri.activate()
        import rpy2.robjects.numpy2ri as numpy2ri
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()

        robjects.r.assign('Mdata',numpy2ri.numpy2ri(numpy_matrix))

        robjects.r('''
            library(Cairo)
            library(ggplot2)
            library(gplots)
            library(RColorBrewer)
            library(vegan)

            Mdata <- Mdata[rowSums(Mdata>0)>0,]
            Mdata <- as.matrix(Mdata)

            core_groups <- length(Mdata[rowSums(Mdata>0)==length(Mdata[1,]),1])
            core_nminus1 <- length(Mdata[rowSums(Mdata>0)>=length(Mdata[1,])-1,1])
            pan_groups <- length(Mdata[,1])

              rarefact <- function(matrice, n_permutatations, sp,ylim){
                  plot(sp,ci.lty=2,xlab="# Genomes", ylab="Orthologous groups", main="", col="blue", ci.col="black",ylim=c(0,ylim))
                  m <- c()
                  for (i in 2:length(matrice[1,])){
                      #print(i)
                      count_list <- c()
                      for (perm in 1:n_permutatations){
                          sub_mat <- matrice[,sample(ncol(matrice), i)]
                          conserved <- sub_mat[rowSums(sub_mat>0)==i,]
                          count_list <- c(count_list,length(conserved[,1]))
                      }

                      #print(rep(i,n_permutatations))
                      #lines(rep(i,n_permutatations),count_list,lty=2, type="l", col="black")
                      segments(i, min(count_list),i, max(count_list),col="blue",lty=2)
                      #print(count_list)
                      m <- c(m,quantile(c(min(count_list),max(count_list)),0.5))
                  }


                  points(seq(2,(length(matrice[2,]))),m,type="l", col="red")
                  print(seq(2:(length(matrice[2,]))))
                  print(m)
                  }
              svg('%s',height=6,width=6)
              if (length(Mdata[1,]) < 10){n_random <- 2} else if (length(Mdata[1,]) < 20) {n_random <- 5} else {n_random <- 10}
              sp <- specaccum(t(Mdata), "random", permutations=n_random)
              rarefact(Mdata, n_random,sp,length(Mdata[,1])*1.3)
              #abline(h=core_groups, col='red')
              #abline(h=pan_groups)
              dev.off()

        ''' % (output_path))

        total = robjects.r['pan_groups']
        core_genome = robjects.r['core_groups']
        core_nminus1 = robjects.r['core_nminus1']
        return total[0], core_genome[0], core_nminus1[0]


def pan_genome_barplot(numpy_matrix,
                       header="",
                       xlab="",
                       ylab="",
                       reverse=False,
                       output_path="~/test.svg"):

        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri
        rpy2.robjects.numpy2ri.activate()
        import rpy2.robjects.numpy2ri as numpy2ri
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()

        robjects.r.assign('Mdata',numpy2ri.numpy2ri(numpy_matrix))

        robjects.r('''
            library(Cairo)
            library(ggplot2)
            library(gplots)
            library(RColorBrewer)
            library(vegan)

            Mdata <- Mdata[rowSums(Mdata>0)>0,]
            Mdata <- as.matrix(Mdata)

            row_counts <- rowSums(Mdata>0)
            t <- table(row_counts)
              svg('%s',height=6,width=6)
                barplot(t,xlab="Number of Genomes", ylab="Number of orthologous groups")
              dev.off()

        ''' % (output_path))

        n_genome_counts = robjects.r['t']

        return n_genome_counts