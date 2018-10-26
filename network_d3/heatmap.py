#!/usr/bin/env python


import datetime
import numpy as np
#import plotly.plotly as py
import plotly.graph_objs as go
import plotly.offline as py
from scipy.spatial.distance import pdist, squareform
import scipy.cluster.hierarchy as hc
import plotly.figure_factory as ff

base = datetime.datetime.today()
date_list = [base - datetime.timedelta(days=x) for x in range(0, 180)]

z = []

mat = '/home/tpillone/work/projets/2018_10_MYCOST/2018_07_4_strains/typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp.tsv'
# /home/tpillone/work/projets/2018_10_MYCOST/2018_07_4_strains/typing/freebayes_joint_genotyping/cgMLST/bwa
import numpy
import networkx
import pandas
import json
import matplotlib.pyplot as plt

m = pandas.read_csv(mat, delimiter='\t', header=0, index_col=0)
link = hc.linkage(m.values, method='centroid')
o1 = hc.leaves_list(link)

mat = m.iloc[o1,:]
mat = mat.iloc[:, o1[::-1]]

nodes = ['S'+ str(i) for i in mat.index]

data = ff.create_annotated_heatmap(
        z= mat.values, # squareform(m.values)
        x= nodes,
        y=nodes,
        colorscale='Reds'
    )

'''
layout = go.Layout(
    title='GitHub commits per day',
    xaxis = dict(ticks='', nticks=36),
    yaxis = dict(ticks='' )
)
'''
#fig = go.Figure(data=data, layout=layout)
py.plot(data, filename="heatmap.html")