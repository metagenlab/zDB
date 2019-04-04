#!/usr/bin/python

import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri
pandas2ri.activate()


def aa_composition_pca(numpy_matrix, target_psoition, path):
        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri as numpy2ri
        rpy2.robjects.numpy2ri.activate()


        robjects.r.assign('Mdata',  numpy2ri.numpy2ri(numpy_matrix))

        robjects.r('''
            print(head(Mdata))

             plot_scores<-function(scores,x,y, target){
                 plot(scores[,x],scores[,y],xlab=paste("comp.",as.character(x)),ylab=paste("comp.",as.character(y)), xlim=range(scores[,c(x,y)]),ylim=range(scores[,c(x,y)]),cex=1.5,pch=20)
                     points(scores[target,x],scores[target,y],pch=18, col="red")
                     #text(scores[25,x],scores[25,y],labels="test",col="red",cex=0.9)
                     abline(h=0,col=2)
                     abline(v=0,col=2)
             }

             visual<-function(groups,clustertable, target){
                 pca2 <- princomp(clustertable)

                 par(mfrow=c(1,3),pty="s")
                 plot_scores(pca2$scores,1,2, target)
                 plot_scores(pca2$scores,1,3,target)
                 plot_scores(pca2$scores,2,3, target)
             }
             png("%s", height=500, width=1300)
             visual(c(1,2,3), Mdata, %s)
             dev.off()
        ''' % (path, target_psoition))

def multiple_aa_composition_pca(numpy_matrix, path):

        '''
        # pca of multiple datasets
        # first column = color factor

        :param numpy_matrix:
        :param target_psoition:
        :param path:
        :return:
        '''

        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri as numpy2ri
        rpy2.robjects.numpy2ri.activate()


        robjects.r.assign('Mdata',  numpy2ri.numpy2ri(numpy_matrix))

        robjects.r('''
                library("FactoMineR")
                library("factoextra")
                print(class(Mdata))
                Mdata[is.na(Mdata)] <- 0
                mat <- as.data.frame(Mdata[,2:length(Mdata[1,])])
                print(head(data.matrix(mat)))
                aa.pca <- PCA(data.matrix(mat), graph = FALSE)
                png("%s", height=600, width=600)
                print(fviz_pca_ind(aa.pca,  label="none", habillage=as.factor(Mdata[,1]))) # ,  label="none", habillage=Mdata[,1]
                dev.off()
        ''' % (path))

import numpy
from biosqldb import manipulate_biosqldb


