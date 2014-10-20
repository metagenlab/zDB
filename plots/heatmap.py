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

def _begin(output,format,new,ratio=1.375,**kwargs):
    """Initializes the plot in *R*."""
    if new:
        if output is None:
            output = unique_filename_in()
        if format == 'pdf':
            robjects.r('pdf("%s",paper="a4",height=8*%f,width=8)' %(output,ratio))
        elif format == 'png':
            robjects.r('png("%s",height=800*%f,width=800,type="cairo")' %(output,ratio))
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
            M[x, y] = M[x, y] + random.uniform(0.01, -0.01)
    return M

def heatmap(M,output=None,format='pdf',new=True,breaks=None, last=True,
            rows=None,columns=None,orderRows=True,orderCols=True,
            return_rowInd=False,cor=False,**kwargs):
    
    """Creates a heatmap of the matrix *M* using *rows* as row labels and *columns* as column labels.
    If either *orderRows* or *orderCols* is True, will cluster accordingly and display a dendrogram."""
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    robjects.r.assign('Mdata',numpy2ri.numpy2ri(M))
    if rows is not None:
        robjects.r.assign('labRow',numpy2ri.numpy2ri(rows))
        plotopt += ",labRow=labRow"
    if columns is not None:
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

myColors=c("white","blue","red")
#myColors = rev(colorRampPalette(brewer.pal(%i,"RdYlBu"))(length(myBreaks)-1))    
par(cex.main=1,oma=c(0,5,0,15))
heatmap.2(as.matrix(Mdata),
          col=myColors, trace="none", breaks=myBreaks, key="False",
          na.rm=TRUE, density.info='none'%s)

legend("bottom", y=NULL, legend=c("single copy", "multicopy"), lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"))



    """ %(ncol, plotopt)) # 
    _end("",last,**kwargs)
    if return_rowInd:
        return (output,array(robjects.r("odhrow")))
    else:
        return output
