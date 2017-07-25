"""

Copy from: https://github.com/bbcf/bbcflib/blob/master/bbcflib/gfminer/figure/rplots.py

These functions use `rpy2` to bind *R* plotting functions with data from numpy arrays.
Each function takes the following arguments:

* output: the filename, a random name will be generated if this is None (default None),
* format: the image format, 'pdf' (default) or 'png',
* new: boolean indicating if a new figure must be started (default) or if the plot is added to the current figure,
* last: boolean, if true this is the last plot on this figure, the file will be finalized and closed on return.
* **kwargs: additional parameters for *R*, such as 'xlab', 'ylab', 'main, 'mfrow', 'log', 'legend'.
"""

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
#from itertools import combinations
from numpy import asarray, array, ndarray, append
import rpy2.robjects.numpy2ri
import random
rpy2.robjects.numpy2ri.activate()

def _begin(output,format,new,height=7,width=14,**kwargs):
    """Initializes the plot in *R*."""
    if new:
        if format == 'pdf':
                robjects.r('library(Cairo);CairoPDF("%s",paper="a4",height=%s,width=%s)' %(output, height, width))
        elif format == 'png':
            robjects.r('png("%s",height=%s*100,width=%s*100,type="cairo")' %(output, height, width))
        elif format == 'svg':
            robjects.r('svg("%s",height=%s*1.2,width=%s*1.2)' %(output, height, width))
        else:
            raise ValueError("Format not supported: %s" %format)
        pars = "lwd=2,cex=1.1,cex.main=1.5,cex.lab=1.3,cex.axis=1.1,mar=c(4,4,1,1),las=1,pch=20"
        if len(kwargs.get('mfrow',[])) == 2:
            pars += ",mfrow=c(%i,%i)" %tuple(kwargs['mfrow'])
        robjects.r('par(%s)' %pars)
    opts = ''
    if 'log' in kwargs: opts += ',log="%s"' %kwargs['log']
    if 'xlim' in kwargs: opts += ',xlim=c(%f,%f)' %tuple(kwargs['xlim'])
    if 'ylim' in kwargs: opts += ',ylim=c(%f,%f)' %tuple(kwargs['ylim'])
    opts += ',main="%s"' %kwargs.get('main','')
    opts += ',xlab="%s"' %kwargs.get('xlab','')
    opts += ',ylab="%s"' %kwargs.get('ylab','')
    return opts, output




def _end(lopts,last,**kwargs):
    """Adds the legend and closes the figure."""
    if not(last): return
    if 'legend' in kwargs:
        names = kwargs.get('legend')
        if names:
            robjects.r("legend(x='topright', legend=%s, col=1:%i%s)" %(list2r(names),len(names),lopts))
    robjects.r("dev.off()")



def randomize_table(M):

    for x in range(0, len(M[:,0])):
        for y in range(0, len(M[0,:])):
            value = float(M[x, y]) + random.uniform(0.01, -0.01)
            M[x, y] = value
    return M




def heatmap_refgenome(M,output=None,format='pdf',new=True,breaks=None, last=True,
            rows=None,columns=None,return_rowInd=False,col_index=[],**kwargs):

    """ create the equivalent of a flat circos:
    on genome as a reference, and ordered columns of presence/absence of homologs"""

    w = len(M[1,:])/float(150) + 16
    h = len(M[:,1]) + 4
    w = h*1.5
    cex = h*0.2
    oma_right = w*2

    plotopt, output = _begin(output=output,format=format,new=new,width=w,height=h,**kwargs)

    robjects.r.assign('Mdata', numpy2ri.numpy2ri(M))

    if len(col_index)>0:
        robjects.r.assign('col_index', col_index)

    if rows is not None:
        import numpy as np
        rows = np.asarray(rows, dtype='a50')#np.chararray(rows)
        robjects.r.assign('labRow',numpy2ri.numpy2ri(rows))
        plotopt += ",labRow=labRow"
    if columns is not None:
        columns = np.asarray(columns, dtype='a50')
        robjects.r.assign('labCol', numpy2ri.numpy2ri(columns))

    if breaks is not None:
        robjects.r("myBreaks = c(%s)" % breaks)
    else:
        robjects.r("myBreaks = seq(ymin,ymax,length.out=15)")
    if kwargs.get('ymax') is not None:
        robjects.r("ymax=%f" %float(kwargs['ymax']))
    else:
        robjects.r("ymax=ceiling(max(Mdata,na.rm=T))")
    try:
        ncol = int(kwargs.get('nb_colors'))
    except (ValueError, TypeError):
        ncol = 10
    ncol = max(3,ncol)
    robjects.r("""

library(gplots)
library(RColorBrewer)
print(myBreaks)
library(vegan)

Mdata[Mdata>1] <- 1

Mdata[,unlist(col_index)] <- Mdata[,unlist(col_index)]*2


#w <- which(colSums(Mdata) != 1)
#Mdata <- Mdata[,w]

print('Drawing heatmap')

myColors=c("white","black","red")
#myColors = rev(colorRampPalette(brewer.pal(%i,"RdYlBu"))(length(myBreaks)-1))
par(oma = c(2, 2, 2, 2), xpd=TRUE)
par(mar = c(1,1,1,1))
par(cex.main=1,oma=c(24,2,2,%s), xpd=TRUE, new=TRUE)

heatmap.2(as.matrix(Mdata),col=myColors, trace="none", breaks=myBreaks, key="False",
 na.rm=TRUE, density.info='none', Rowv=TRUE, Colv=FALSE,labRow=labRow,labCol=labCol,cexRow=6,cexCol=4
 ,distfun=function(m) vegdist(m,method="jaccard", binary=TRUE)) # , Colv=cdendo

#legend("topright", inset=c(0.4,1.5), y=NULL, legend=c("single copy", "multicopy"), lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"))



    """ % (ncol, oma_right)) #
    _end("",last,**kwargs)

    if return_rowInd:
        return (output,array(robjects.r("odhrow")))
    else:
        return output


def heatmap_pangenome(M,output=None,format='pdf',new=True,breaks=None, last=True,
            rows=None,columns=None,orderRows=True,orderCols=True,
            return_rowInd=False,cor=False,**kwargs):

    """Creates a heatmap of the matrix *M* using *rows* as row labels and *columns* as column labels.
    If either *orderRows* or *orderCols* is True, will cluster accordingly and display a dendrogram."""

    w = len(M[1,:])/float(150) + 16
    h = len(M[:,1]) + 6
    w = h*1.5
    cex = h*0.2
    oma_right = w*3

    print 'size', h, w, oma_right

    plotopt, output = _begin(output=output,format=format,new=new,width=w,height=h,**kwargs)
    print "ok"
    robjects.r.assign('Mdata', numpy2ri.numpy2ri(M))
    if rows is not None:
        import numpy as np
        rows = np.asarray(rows, dtype='a50')#np.chararray(rows)
        robjects.r.assign('labRow',numpy2ri.numpy2ri(rows))
        plotopt += ",labRow=labRow"
    if columns is not None:
        columns = np.asarray(columns, dtype='a50')
        robjects.r.assign('labCol', numpy2ri.numpy2ri(columns))
        plotopt += ",labCol=labCol"

    if return_rowInd:
        orderRows = True
        robjects.r("""
if (exists("myCor")){hrow = as.dendrogram(hclust(myCor(Mdata)))} else {hrow = as.dendrogram(hclust(dist(Mdata)))}
odhrow = rev(order.dendrogram(hrow))
labRow = rep('',nrow(Mdata))
nb_IDs = ceiling(nrow(Mdata)/30)
labRow[seq(1,nrow(Mdata),nb_IDs)] = seq(1,nrow(Mdata),nb_IDs)
labRow[odhrow] = labRow""")
        plotopt += ",Rowv=hrow"
    if orderCols and orderRows:
        plotopt += ",dendrogram='both',lhei=c(2,10,1,2),lwid=c(1,3),mar=c(2,2),lmat=matrix(c(0,2,0,0,3,1,0,4),ncol=2)"
    elif orderCols:
        plotopt += ",Rowv=F,dendrogram='column',lhei=c(2,10,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(3,1,2,4),ncol=1)"
    elif orderRows:
        plotopt += ",Colv=F,dendrogram='row',lhei=c(10,1,2),lwid=c(1,3),mar=c(2,2),lmat=matrix(c(2,0,3,1,0,4),ncol=2)"
    else:
        plotopt += ",Colv=F,Rowv=F,dendrogram='none',lhei=c(10,1,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(1,2,3,4),ncol=1)"
    plotopt += ',cexRow=%s, cexCol=%s' % (cex, cex)
    if kwargs.get('ymin') is not None:
        robjects.r("ymin=%f" %float(kwargs['ymin']))
    else:
        robjects.r("ymin=floor(min(Mdata,na.rm=T))")
    if breaks is not None:
        robjects.r("myBreaks = c(%s)" % breaks)
    else:
        robjects.r("myBreaks = seq(ymin,ymax,length.out=15)")
    if kwargs.get('ymax') is not None:
        robjects.r("ymax=%f" %float(kwargs['ymax']))
    else:
        robjects.r("ymax=ceiling(max(Mdata,na.rm=T))")
    try:
        ncol = int(kwargs.get('nb_colors'))
    except (ValueError, TypeError):
        ncol = 10
    ncol = max(3,ncol)
    robjects.r("""

library(gplots)
library(RColorBrewer)
print(myBreaks)
library(vegan)
Mdata[Mdata>1] <- 1

#w <- which(colSums(Mdata) != 1)
#Mdata <- Mdata[,w]

print('Drawing heatmap')

myColors=c("white","darkgray","black")
#myColors = rev(colorRampPalette(brewer.pal(%i,"RdYlBu"))(length(myBreaks)-1))
par(oma = c(12, 5, 5, 5), xpd=TRUE)
par(mar = c(1,1,1,1))
par(cex.main=1,oma=c(15,5,5,%s), xpd=TRUE, new=TRUE)

cd <- vegdist(as.matrix(t(Mdata)), method = "jaccard", binary=TRUE)
cfit <- hclust(cd, method="average") # , method="average"
Colv <- colSums(as.matrix(Mdata), na.rm = TRUE) # rep(0,length(Mdata[,1]))
cdendo <- reorder(as.dendrogram(cfit), Colv)

m <- match(labels(cdendo), seq(1:length(Mdata[1,])))

if (length(cfit$order) < 3000) {

    labels <- rep("",length(cfit$order))
    for (i in seq(0,length(cfit$order),100)){
        labels[m[i]] <- i
        }
    labels[m[length(cfit$order)]] <- length(cfit$order)

}else{

    labels <- rep("",length(cfit$order))
    for (i in seq(0,length(cfit$order),1000)){
        labels[m[i]] <- i
        }
    labels[m[length(cfit$order)]] <- length(cfit$order)
}


heatmap.2(as.matrix(Mdata),col=myColors, trace="none", breaks=myBreaks, key="False",distfun=function(m) vegdist(m,method="jaccard", binary=TRUE),
          na.rm=TRUE, density.info='none'%s, labCol=labels, hclustfun=function(m) hclust(m, method="average")) # , Colv=cdendo

#legend("topright", inset=c(0.4,1.5), y=NULL, legend=c("single copy", "multicopy"), lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"))



    """ %(ncol, oma_right, plotopt)) #
    _end("",last,**kwargs)

    if return_rowInd:
        return (output,array(robjects.r("odhrow")))
    else:
        return output





def heatmap(M,output=None,format='pdf',new=True,breaks=None, last=True,
            rows=None,columns=None,orderRows=True,orderCols=True,
            return_rowInd=False,cor=False,**kwargs):
    
    """Creates a heatmap of the matrix *M* using *rows* as row labels and *columns* as column labels.
    If either *orderRows* or *orderCols* is True, will cluster accordingly and display a dendrogram."""

    print M
    print "cols", columns
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    print "ok"
    robjects.r.assign('Mdata',numpy2ri.numpy2ri(M))
    if rows is not None:
        import numpy as np
        rows = np.chararray(rows)
        robjects.r.assign('labRow',numpy2ri.numpy2ri(rows))
        plotopt += ",labRow=labRow"
    if columns is not None:
        import numpy as np
        print "cols", columns
        columns = np.asarray(columns, dtype='a50')
        print (len(columns))
        robjects.r.assign('labCol',numpy2ri.numpy2ri(columns))
        plotopt += ",labCol=labCol"
    if cor:

        robjects.r("myCor = function(x) {as.dist(1-cor(t(x),use='pairwise.complete.obs'))}")
        plotopt += ", distfun=myCor"
    if return_rowInd:
        orderRows = True
        robjects.r("""
if (exists("myCor")){hrow = as.dendrogram(hclust(myCor(Mdata)))} else {hrow = as.dendrogram(hclust(dist(Mdata)))}
odhrow = rev(order.dendrogram(hrow))
labRow = rep('',nrow(Mdata))
nb_IDs = ceiling(nrow(Mdata)/30)
labRow[seq(1,nrow(Mdata),nb_IDs)] = seq(1,nrow(Mdata),nb_IDs)
labRow[odhrow] = labRow""")
        plotopt += ",Rowv=hrow"
    if orderCols and orderRows:
        plotopt += ",dendrogram='both',lhei=c(2,10,1,2),lwid=c(1,3),mar=c(2,2),lmat=matrix(c(0,2,0,0,3,1,0,4),ncol=2)"
    elif orderCols:
        plotopt += ",Rowv=F,dendrogram='column',lhei=c(2,10,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(3,1,2,4),ncol=1)"
    elif orderRows:
        plotopt += ",Colv=F,dendrogram='row',lhei=c(10,1,2),lwid=c(1,3),mar=c(2,2),lmat=matrix(c(2,0,3,1,0,4),ncol=2)"
    else:
        plotopt += ",Colv=F,Rowv=F,dendrogram='none',lhei=c(10,1,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(1,2,3,4),ncol=1)"
    if kwargs.get('ymin') is not None:
        robjects.r("ymin=%f" %float(kwargs['ymin']))
    else:
        robjects.r("ymin=floor(min(Mdata,na.rm=T))")
    if breaks is not None:
        robjects.r("myBreaks = c(%s)" % breaks)
    else:
        robjects.r("myBreaks = seq(ymin,ymax,length.out=15)")
    if kwargs.get('ymax') is not None:
        robjects.r("ymax=%f" %float(kwargs['ymax']))
    else:
        robjects.r("ymax=ceiling(max(Mdata,na.rm=T))")
    try:
        ncol = int(kwargs.get('nb_colors'))
    except (ValueError, TypeError):
        ncol = 10
    ncol = max(3,ncol)
    robjects.r("""

library(gplots)
library(RColorBrewer)
print(myBreaks)
print('Drawing heatmap')
myColors=c("white","blue","red")
#myColors = rev(colorRampPalette(brewer.pal(%i,"RdYlBu"))(length(myBreaks)-1))
par(oma = c(22, 0, 0, 5), xpd=TRUE)
par(mar = c(1,1,1,1))
par(cex.main=1,oma=c(22,0,0,5), xpd=TRUE, new=TRUE)
heatmap.2(as.matrix(Mdata),
          col=myColors, trace="none", breaks=myBreaks, key="False",
          na.rm=TRUE, density.info='none'%s)

#legend("topright", inset=c(0.4,1.5), y=NULL, legend=c("single copy", "multicopy"), lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"))



    """ %(ncol, plotopt)) #
    _end("",last,**kwargs)
    if return_rowInd:
        return (output,array(robjects.r("odhrow")))
    else:
        return output

def heatmap_ksnp(M,output=None,format='pdf',new=True,breaks=None, last=True,
            rows=None,columns=None,orderRows=True,orderCols=True,
            return_rowInd=False,cor=False,**kwargs):

    """Creates a heatmap of the matrix *M* using *rows* as row labels and *columns* as column labels.
    If either *orderRows* or *orderCols* is True, will cluster accordingly and display a dendrogram."""
        
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    robjects.r.assign('Mdata',numpy2ri.numpy2ri(M))
    robjects.r("""
    rowsums <- rowSums(Mdata)
    w <-which(rowsums<50)
    Mdata <-Mdata[w,]


    """)
    if rows is not None:
        robjects.r.assign('labRow',numpy2ri.numpy2ri(rows))
        plotopt += ",labRow=labRow"
    if columns is not None:
        print "columns", columns
        import numpy as np
        columns = np.asarray(columns, dtype='a50')
        robjects.r.assign('labCol',numpy2ri.numpy2ri(columns))
        plotopt += ",labCol=labCol"
    if cor:
        robjects.r("myCor = function(x) {as.dist(1-cor(t(x),use='pairwise.complete.obs'))}")
        plotopt += ", distfun=myCor"
    if return_rowInd:
        orderRows = True
        robjects.r("""
if (exists("myCor")){hrow = as.dendrogram(hclust(myCor(Mdata)))} else {hrow = as.dendrogram(hclust(dist(Mdata)))}
odhrow = rev(order.dendrogram(hrow))
labRow = rep('',nrow(Mdata))
nb_IDs = ceiling(nrow(Mdata)/30)
labRow[seq(1,nrow(Mdata),nb_IDs)] = seq(1,nrow(Mdata),nb_IDs)
labRow[odhrow] = labRow""")
        plotopt += ",Rowv=hrow"
    if orderCols and orderRows:
        plotopt += ",dendrogram='both',lhei=c(2,10,1,2),lwid=c(1,3),mar=c(2,2),lmat=matrix(c(0,2,0,0,3,1,0,4),ncol=2)"
    elif orderCols:
        plotopt += ",Rowv=F,dendrogram='column',lhei=c(2,10,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(3,1,2,4),ncol=1)"
    elif orderRows:
        plotopt += ",Colv=F,dendrogram='row',lhei=c(10,1,2),lwid=c(1,3),mar=c(2,2),lmat=matrix(c(2,0,3,1,0,4),ncol=2)"
    else:
        plotopt += ",Colv=F,Rowv=F,dendrogram='none',lhei=c(10,1,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(1,2,3,4),ncol=1)"
    if kwargs.get('ymin') is not None:
        robjects.r("ymin=%f" %float(kwargs['ymin']))
    else:
        robjects.r("ymin=floor(min(Mdata,na.rm=T))")
    if breaks is not None:
        robjects.r("myBreaks = c(%s)" % breaks)
    else:
        robjects.r("myBreaks = seq(ymin,ymax,length.out=15)")
    if kwargs.get('ymax') is not None:
        robjects.r("ymax=%f" %float(kwargs['ymax']))
    else:
        robjects.r("ymax=ceiling(max(Mdata,na.rm=T))")
    try:
        ncol = int(kwargs.get('nb_colors'))
    except (ValueError, TypeError):
        ncol = 10
    ncol = max(3,ncol)
    robjects.r("""
library(gplots)
library(RColorBrewer)
print(myBreaks)
myColors=c("blue","red","white")
#myColors = rev(colorRampPalette(brewer.pal(%i,"RdYlBu"))(length(myBreaks)-1))
par(oma = c(22, 0, 0, 5), xpd=TRUE)
par(mar = c(1,1,1,1))
par(cex.main=1,oma=c(22,0,0,5), xpd=TRUE, new=TRUE)
heatmap.2(as.matrix(Mdata),
          col=myColors, trace="none", breaks=myBreaks, key="False",
          na.rm=TRUE, density.info='none'%s)

#legend("topright", inset=c(0.4,1.5), y=NULL, legend=c("single copy", "multicopy"), lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"))



    """ %(ncol, plotopt)) #
    _end("",last,**kwargs)
    if return_rowInd:
        return (output,array(robjects.r("odhrow")))
    else:
        return output
