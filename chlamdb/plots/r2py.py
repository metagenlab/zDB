#!/usr/bin/python

from numpy import asarray, array, ndarray, append
import rpy2.robjects.numpy2ri
import random
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
rpy2.robjects.numpy2ri.activate()
from manipulate_trees import biodb2heatmap

a = array([[1,2,2,2,2,2,2],[1,1,1,1,1,1,1],[1,2,1,2,1,2,1]])

robjects.r.assign('Mdata',numpy2ri.numpy2ri(a))
output = "tata.pdf"
ratio = 0.5

robjects.r('pdf("%s",paper="a4",height=8*%f,width=8)' %(output,ratio))
robjects.r('plot(Mdata)')
robjects.r('dev.off()')

biodb2heatmap("chlamydia_02_15")