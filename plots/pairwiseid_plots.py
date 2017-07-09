#!/usr/bin/python

import rpy2.robjects.numpy2ri

import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr
import rpy2.rlike.container as rlc
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import pandas.rpy.common as com



def density_plot(value_list_of_lists, label_list,
                 header="-",
                 xlab="-",
                 ylab="-",
                 output_path="~/test.svg",
                 type="hist",
                 abline_list=[],
                 show_median=True,
                 min_value=False,
                 max_value=False,
                 breaks_manual=False):

        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri
        rpy2.robjects.numpy2ri.activate()
        import rpy2.robjects.numpy2ri as numpy2ri
        from pandas import DataFrame
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()


        format_list = []
        for value_list, label in zip(value_list_of_lists, label_list):
            for value in value_list:
                format_list.append([value, label])
        df = DataFrame(format_list, columns=['identity', 'comp'])
        #m = m.astype(float)
        robjects.r.assign('plot_data', pandas2ri.py2ri(df))

        if show_median:
            median = ' + geom_vline(data=mu, aes(xintercept=identity.mean, color=comp), linetype="dashed")'

        else:
            median = ""

        if min_value:
            limits = '+ scale_x_continuous(limits = c(%s, %s),breaks = seq(%s, %s, 5))' % (min_value,
                                                                                            max_value,
                                                                                            min_value,
                                                                                            max_value)
        else:
            limits = ''

        if breaks_manual:
            breaks = 'breaks=pretty(range(plot_data$identity), n = %s, min.n = 1),' % breaks_manual
        else:
            breaks= ''

        if type == 'density':
            plot_code = '''

           d <- density(as.numeric(value_list))
           w <- which(d$y == max(d$y))
           plot(d, main="%s", xlab="%s", ylab="%s", xlim=c(0,100))
           abline(v=d$x[w], col="red")
            '''% (
               header,
               xlab,
               ylab)
        elif type == 'hist2':
            plot_code = '''

           hist(as.numeric(value_list), main="%s", xlab="%s", ylab="%s", xlim=c(0,100))

            '''% (
               header,
               xlab,
               ylab)
        elif type == "hist":
            plot_code = '''

            p <- ggplot(plot_data, aes(identity, fill = comp)) + geom_histogram(%s alpha = 0.5, aes(y = ..density..), position = 'identity')
            p <- p + geom_density(alpha=0.01, aes(identity, colour = comp)) %s
            p <- p + xlab("%s") +  ylab("%s") + ggtitle("%s")+theme(legend.position="bottom")

            ''' % (breaks,median,
               xlab,
               ylab,
                   header)
        else:
            print "unkown type"

        robjects.r.assign('abline_list', abline_list)

        if len(abline_list)>0:
            abline = "+ geom_vline(xintercept = unlist(abline_list))"
        else:
            abline = ""

        robjects.r.assign('value_list_of_lists', value_list_of_lists)
        robjects.r.assign('label_list', label_list)

        robjects.r('''
            library(genoPlotR)
            library(Cairo)
            library(ggplot2)

            library(plyr)
            mu <- ddply(plot_data, "comp", summarise, identity.mean=median(identity))

            plot_data$identity <- as.numeric(plot_data$identity)

            CairoSVG('%s',height=8,width=8)
            p <- %s
            print (p + theme_bw() %s +theme(axis.text=element_text(size=10),
            axis.title=element_text(size=14))+ theme(legend.text=element_text(size=10)) %s
             +theme(legend.position="bottom",legend.direction="vertical"))
            dev.off()


        ''' % (output_path, plot_code, abline, limits))

def basic_plot(values_x, values_y=False, header="", xlab="", ylab="", output_path="~/test.svg", type="hist"):

        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri
        rpy2.robjects.numpy2ri.activate()
        import rpy2.robjects.numpy2ri as numpy2ri
        from pandas import DataFrame
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()

        a = np.asarray(values_x, dtype='float')
        if values_y:
            b = np.asarray(values_y, dtype='float')

            robjects.r.assign('values_x', numpy2ri.numpy2ri(a))
            robjects.r.assign('values_y', numpy2ri.numpy2ri(b))

            robjects.r('''
                library(genoPlotR)
                library(Cairo)
                library(ggplot2)


                library(plyr)
                mu <- ddply(plot_data, "comp", summarise, identity.mean=median(identity))
                print (mu)
                plot_data$identity <- as.numeric(plot_data$identity)

                svg('%s',height=6,width=14)
                plot(values_x, values_y, pch=20) # , ylim=c(0,100)
                dev.off()


            ''' % (output_path))
        else:
            robjects.r.assign('values_x', numpy2ri.numpy2ri(a))

            robjects.r('''
                library(genoPlotR)
                library(Cairo)
                library(ggplot2)

                svg('%s',height=7,width=7)
                barplot(table(values_x), main="Conservation of predicted effectors in other genomes")
                dev.off()

            ''' % (output_path))

def identity_heatmap_plot(numpy_matrix, labels,
                          header="",
                          xlab="",
                          ylab="",
                          reverse=False,
                          output_path="~/test.svg"):

        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri
        rpy2.robjects.numpy2ri.activate()
        import rpy2.robjects.numpy2ri as numpy2ri
        from pandas import DataFrame
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()

        robjects.r.assign('Mdata',numpy2ri.numpy2ri(numpy_matrix))
        b = np.asarray(labels, dtype='str')
        robjects.r.assign('labels',numpy2ri.numpy2ri(b))

        if reverse:
            plot = '''
                cols <- rev(brewer.pal(9, "Blues"))
                heatmap.2(100-as.matrix(Mdata), trace="none", key="True",col=cols,
                na.rm=TRUE, density.info='none', cellnote=as.matrix(Mdata), notecol="black", labRow=labels, labCol=labels)
            '''
        else:
            plot = '''
                cols <- brewer.pal(9, "Blues")
                heatmap.2(as.matrix(Mdata), trace="none", key="True",col=cols,
                na.rm=TRUE, density.info='none', cellnote=as.matrix(Mdata), notecol="black", labRow=labels, labCol=labels)
            '''


        robjects.r('''
            library(Cairo)
            library(ggplot2)
            library(gplots)
            library(RColorBrewer)


            h <- length(Mdata[,1])/4+8
            w <- length(Mdata[,1])/2+8
            print(h)
            svg('%s',height=h,width=w)
                par(oma = c(5, 0, 0, 8), xpd=TRUE)
                par(mar = c(5,1,1,8))
                par(cex.main=1,oma=c(22,0,0,20), xpd=TRUE, new=TRUE)
                %s
            dev.off()

        ''' % (output_path, plot))

