#!/usr/bin/python

import rpy2.robjects.numpy2ri

import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr
import rpy2.rlike.container as rlc
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
#import pandas.rpy.common as com



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
                 breaks_manual=False,
                 show_mode=False):

        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri
        rpy2.robjects.numpy2ri.activate()
        import rpy2.robjects.numpy2ri as numpy2ri
        from pandas import DataFrame
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()

        header_split = header.split(" ")
        if len(header) > 2 and not "bacterium" in header and not "sp." in header:
            header1 = ' '.join(header_split[0:2])
            header2 = ' '.join(header_split[2:])
        else:
            header1 = ' '.join(header_split[0:1])
            header2 = ' '.join(header_split[1:])

        format_list = []
        for value_list, label in zip(value_list_of_lists, label_list):
            for value in value_list:
                format_list.append([value, label])
        df = DataFrame(format_list, columns=['identity', 'comp'])
        #m = m.astype(float)
        robjects.r.assign('plot_data', pandas2ri.py2ri(df))

        median = ""
        if show_median:
            median += ' + geom_vline(data=mu, aes(xintercept=identity.mean, color=comp), linetype="dashed")'

        if show_mode:
            median += ' + geom_vline(data=mu2, aes(xintercept=identity.mode, color=comp), size = 0.5)'


        if min_value is not False:
            limits = '+ scale_x_continuous(limits = c(%s, %s),breaks = seq(%s, %s, 5))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))' % (min_value,
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
            p <- p + xlab("%s") +  ylab("%s") + ggtitle(~ paste(italic("%s"), "%s")) + theme(legend.position="bottom")+ theme_bw()

            ''' % (breaks,
                   median,
                   xlab,
                   ylab,
                   header1, header2)
        else:
            print ("unkown type")

        robjects.r.assign('abline_list', abline_list)

        if len(abline_list)>0:
            abline = "+ geom_vline(xintercept = unlist(abline_list))"
        else:
            abline = ""

        robjects.r.assign('value_list_of_lists', value_list_of_lists)
        robjects.r.assign('label_list', label_list)

        robjects.r('''
            #library(genoPlotR)
            #library(Cairo)
            library(ggplot2)
            #library(flextable)
            library(svglite)


            my_svg <- function(file, width, height) {
              library(RSvgDevice)
              devSVG(file = file, width = width, height = height, onefile = TRUE, xmlHeader = TRUE)
            }


            library(plyr)
            mu <- ddply(plot_data, "comp", summarise, identity.mean=median(identity))
            print(mu)
            plot_data$identity <- as.numeric(plot_data$identity)

            modef <- function(col){

           d <- density(as.numeric(col))
           w <- which(d$y == max(d$y))
            return (d$x[w])
            }

            mu2 <- ddply(plot_data, "comp", summarise, identity.mode=modef(identity))

            

            p <- %s
            p <- p + theme_bw() %s +theme(axis.text=element_text(size=10), axis.title=element_text(size=16)) 
            p <- p + theme(legend.text=element_text(size=13)) %s 
            p <- p + theme(legend.position="bottom",legend.direction="vertical", plot.margin = margin(1, 1, 1, 1, "cm"),)
            p <- p + theme(legend.title=element_blank())+ theme(plot.title = element_text(size=18, face="bold.italic"))
            p <- p #+ scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.01))
            #p <- p + labs(title = ~ paste("test", italic("Species")))
            
            diff <- abs(5-length(label_list))
            h <- 8 - (0.24*diff)
            
            svglite('%s',height=h,width=8)
            print(p)
            dev.off()


        ''' % (plot_code, abline, limits, output_path))

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
                #library(genoPlotR)
                library(Cairo)
                library(ggplot2)


                library(plyr)
                #mu <- ddply(plot_data, "comp", summarise, identity.mean=median(identity))
                #print (mu)
                #plot_data$identity <- as.numeric(plot_data$identity)

                svg('%s',height=6,width=14)
                plot(values_x, values_y, pch=20) # , ylim=c(0,100)
                dev.off()


            ''' % (output_path))
        else:
            robjects.r.assign('values_x', numpy2ri.numpy2ri(a))

            robjects.r('''
                #library(genoPlotR)
                library(Cairo)
                library(ggplot2)

                svg('%s',height=7,width=7)
                #barplot(table(values_x), main="Conservation of predicted effectors in other genomes")

                  library(ggplot2)

                  mytable <- as.data.frame(table(values_x))
                  #print(mytable)
                  p <- ggplot(mytable, aes(x = reorder(values_x, -order(values_x)), y = Freq)) + geom_bar(stat = "identity")
                  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))+ coord_flip()
                  print(p)
                dev.off()

            ''' % (output_path))

def plot_multiseries_points(data,
                            output_path,
                            x_label = "Reciprocal BBH median identity (%)",
                            y_label = "16S rRNA identity (%)",
                            x_label2 = "Identity (%)" ):

        import pandas as pd
        df = pd.DataFrame(data, columns=['category','val_x', 'val_y'])
        #m = m.astype(float)
        robjects.r.assign('plot_data', pandas2ri.py2ri(df))

        robjects.r('''
            #library(genoPlotR)
            library(Cairo)
            library(ggplot2)
            library(svglite)

            print(head(plot_data))

            svglite('%s.svg',height=7,width=9)
            p <- ggplot(data = plot_data, aes(x = val_x, y = val_y, color = category)) + geom_point(alpha = 0.9, shape=2, size=1.2) # alpha = 1, shape=1, size=2
            #p <- p + scale_x_continuous(limits = c(35, 100),breaks = seq(35, 100, 5))
            p <- p + scale_x_continuous(limits = c(30, 100),breaks = seq(30, 100, 5)) # c(0, 1),breaks = seq(0, 1, 0.1)
            p <- p + scale_y_continuous(limits = c(82, 100),breaks = seq(82, 100, 2))
            #p <- p + geom_smooth(method=lm, se=FALSE)
            p <- p + theme_bw() + theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=16), legend.text = element_text(size = 15), legend.key.size = unit(1.5,"line"), plot.margin = margin(1, 1, 1, 1, "cm"))
            p <-  p + guides(color = guide_legend(override.aes = list(size=4)))
            p <- p + xlab("%s") + ylab("%s")
            print(p)
            dev.off()

            svglite('%s_density.svg',height=7,width=8)
            p <- ggplot(plot_data, aes(val_x, fill = category)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
            p <- p + geom_density(alpha=0.01, aes(val_x, colour = category))+ theme_bw()
            p <- p + theme(legend.title=element_blank(),axis.text=element_text(size=10), axis.title=element_text(size=16), legend.text = element_text(size = 15), plot.margin = margin(1, 1, 1, 1, "cm"))
            p <- p + xlab("%s") + ylab("Density")
            print(p)
            dev.off()


        ''' % (output_path, x_label, y_label, output_path, x_label2))





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
                h <- heatmap.2(as.matrix(Mdata), trace="none", key="True",col=cols,
                na.rm=TRUE, density.info='none', cellnote=as.matrix(Mdata), notecol="black", labRow=labels, labCol=labels)
            '''


        robjects.r('''
            library(Cairo)
            library(ggplot2)
            library(gplots)
            library(RColorBrewer)

            rownames(Mdata) <- labels
            colnames(Mdata) <- labels


            h <- length(Mdata[,1])/4+8
            w <- length(Mdata[,1])/2+8
            print(h)
            print(w)
            svg('%s',height=h,width=w)
                par(oma = c(5, 0, 0, 8), xpd=TRUE)
                par(mar = c(5,1,1,8))
                par(cex.main=1,oma=c(22,0,0,20), xpd=TRUE, new=TRUE)
                %s
                #print(h)
                write.table(Mdata, "/home/trestan/identity.tab", sep="\t",col.names = NA)
            dev.off()

        ''' % (output_path, plot))

