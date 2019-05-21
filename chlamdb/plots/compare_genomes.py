#!/usr/bin/env python
# python 2.7.5 requires biopython

########### promer2circos ############

# A script for converting show-coords output for sequence matches from the Mummerplot package to circos link-track format.


# Version 1. Adam Taranto, January 2014.
# Contact Adam Taranto, adam.p.taranto@gmail.com

###Import modules
import sys;
import argparse;
import re;

# coords_input=open('data/promerdata.txt','rU').readlines()
import random
import string









def purge(dir, pattern):
    import os
    for f in os.listdir(dir):
    	if re.search(pattern, f):
    		os.remove(os.path.join(dir, f))

class CircosException(Exception):
    pass


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]


class GapOverlapIdentification():
    '''
    identify gaps from nucmer or any of list of alignments
    input: dictionnay of contigs with a corresponding list of start positions (alignments)
    max_merge_disence: maximal distence between gaps allowed for the merging of close gaps

    '''

    def __init__(self):
        self.gap_data = {}


    def get_gap_data_from_start_stop_list(self,contig2start_stop_lists, max_merge_distence=2000, add_coords = False):

        import pandas as pd
        import numpy as np

        for contig in contig2start_stop_lists:
            # order data with pandas dataframe: order based on start of the alignment
            data = pd.DataFrame({'start': contig2start_stop_lists[contig]["start"],
                    'stop': contig2start_stop_lists[contig]["stop"] })
            data.start = data.start.astype(np.int64)
            data.stop = data.stop.astype(np.int64)
            data_sort = data.sort_values(by=["start"])#sort(columns=["start"])

            # the 2 index of alignments to compare
            index_start = 0
            comparison_index = 1
            # iter alignments
            while index_start < len(data_sort['start']):
                # compare an alignment with all alignments located "en amont"
                # if cumulated index outsize the dataframe, comparison is finished
                if index_start+comparison_index >= len(data_sort['start']):
                    break

                # case if nex alignment included within the curent alignment: skip
                if int(data_sort['stop'][index_start+comparison_index])<= int(data_sort['stop'][index_start]):
                    comparison_index+=1
                    continue
                # case if second frament start within first fragment and is longer than the first one
                # jump to the next alignment (skipping small alignments within bigget alignments)
                elif int(data_sort['start'][index_start+comparison_index])<=int(data_sort['stop'][index_start]) and  \
                    int(data_sort['stop'][index_start+comparison_index]) > int(data_sort['stop'][index_start]):
                    index_start+=comparison_index
                    comparison_index=1
                    continue
                # case if next alignment is not "a cheval" with the current one
                # if distence big enough: gap
                elif int(data_sort['start'][index_start+comparison_index])-int(data_sort['stop'][index_start]) > max_merge_distence:
                    # create a distionnary to store gap data
                    if contig not in self.gap_data:
                        self.gap_data[contig] = [[data_sort['stop'][index_start], data_sort['start'][index_start+comparison_index]]]
                    else:
                        self.gap_data[contig].append([data_sort['stop'][index_start], data_sort['start'][index_start+comparison_index]])

                    index_start+=comparison_index
                    comparison_index=1
                else:
                    #print 'add'
                    index_start+=comparison_index
                    comparison_index=1
        return self.gap_data


    def get_circos_file_from_gap_data(self, out_highlight, out_labels):
        import numpy as np

        gap_data = self.gap_data

        h = open(out_highlight, 'w')
        g = open(out_labels, 'w')

        for contig in gap_data:
            print gap_data[contig]
            if len(gap_data[contig])>1:
                start_list = []
                stop_list = []
                for gap in gap_data[contig]:
                    start_list.append(gap[0])
                    stop_list.append(gap[1])
                print "start_list"
                print start_list
                print "stop_list"
                print stop_list
                data = pd.DataFrame({'start': start_list,
                'stop': stop_list })
                data.start = data.start.astype(np.int64)
                data.stop = data.stop.astype(np.int64)
                #print 'add', data['start'][0]+1
                data_sort = data.sort(columns=["start"])
                print data_sort
                start = data_sort['start'][0]
                stop = data_sort['stop'][0]

                for i in range(0, len(data_sort['start'])-1):
                    print i, len(data_sort['start'])
                    # if less than 10kb between 2 gaps, merge
                    if data_sort['start'][i+1]-data_sort['stop'][i] <10000:
                        stop = data_sort['stop'][i+1]
                        # if the last gap considered, write it
                        if (i+2) == len(data_sort['start']):
                            # gap highlight
                            h.write("%s\t%s\t%s\n" % (contig,start, stop))
                            # gap lables (start and stop)
                            g.write("%s\t%s\t%s\t%s\n"% (contig,start, start+5, start))
                            g.write("%s\t%s\t%s\t%s\n"% (contig,stop, stop+5, stop))
                    # more than 10kb between the 2 gaps
                    else:
                        if i+2 == len(data_sort['start']):
                            # last gaps compared. as not merged, write both
                            h.write("%s\t%s\t%s\n" % (contig,start, stop))
                            # write labels
                            g.write("%s\t%s\t%s\t%s\n"% (contig,start, start+5, start))
                            g.write("%s\t%s\t%s\t%s\n"% (contig,stop, stop+5, stop))
                            h.write("%s\t%s\t%s\n" % (contig,data_sort['start'][i+1], data_sort['stop'][i+1]))
                            # write labels
                            g.write("%s\t%s\t%s\t%s\n"% (contig,data_sort['start'][i+1], data_sort['start'][i+1]+5, data_sort['start'][i+1]))
                            g.write("%s\t%s\t%s\t%s\n"% (contig,data_sort['stop'][i+1], data_sort['stop'][i+1]+5, data_sort['stop'][i+1]))
                        else:
                            # not merged, write gap and reinitialize start and stop of the next gap
                            h.write("%s\t%s\t%s\n" % (contig,start, stop))
                            # write labels
                            g.write("%s\t%s\t%s\t%s\n"% (contig,start, start+5, start))
                            g.write("%s\t%s\t%s\t%s\n"% (contig,stop, stop+5, stop))
                            start = data_sort['start'][i+1]
                            stop = data_sort['stop'][i+1]

                #print data_sort
                #import sys
                #sys.exit()


            else:
                start = gap_data[contig][0][0]
                stop = gap_data[contig][0][1]
                h.write("%s\t%s\t%s\n" % (contig, start, stop))
                # write labels
                g.write("%s\t%s\t%s\t%s\n"% (contig,start, start+5, start))
                g.write("%s\t%s\t%s\t%s\n"% (contig,stop, stop+5, stop))


        h.close()
        g.close()



class CompareWholeGenomes():

    def __init__(self):
        pass

    def justLinks(self, coords_input):
        ##Find first row after header
        i = 0
        for row in coords_input:
            if len(re.split(r'\t+', row)) > 3:
                headerRow = i + 1
                break
            else:
                i += 1

        return coords_input[headerRow:None]

    def execute_promer(self,fasta1, fasta2, out_file = "out.coords", algo="nucmer"):
        import shell_command

        #cmd1 = 'promer --mum -l 5 %s %s' % (fasta1, fasta2)
        for one_fasta in fasta2:
            if algo == 'nucmer':
                cmd1 = 'nucmer -mum -b 200 -c 65 -g 90 -l 20 %s %s' % (fasta1, one_fasta)
                cmd2 = 'show-coords -T -r -c -L 100 -I 30 out.delta > %s.coords' % one_fasta.split('.')[0]
                a, b, c = shell_command.shell_command(cmd1)
                a, b, c = shell_command.shell_command(cmd2)
            elif algo == 'promer':
                cmd1 = 'promer --mum -l 5 %s %s' % (fasta1, one_fasta)
                cmd2 = 'show-coords -T -r -c -L 100 -I 30 out.delta > %s.coords' % one_fasta.split('.')[0]
                a, b, c = shell_command.shell_command(cmd1)
                a, b, c = shell_command.shell_command(cmd2)


    def execute_megablast(self,fasta1, fasta2):
        import shell_command
        for one_fasta in fasta2:
            cmd1 = "formatdb -i %s -p F" % one_fasta
            cmd2 = 'blastn -task megablast -query %s -db %s -evalue 1e-5 -outfmt 6 -out blast_result_%s.tab' % (fasta1,
                                                                                                             one_fasta,
                                                                                                             one_fasta.split('.')[0])
            a, b, c = shell_command.shell_command(cmd1)
            a, b, c = shell_command.shell_command(cmd2)


    def nucmer_coords2data(self, coords_input, algo="nucmer", contigs_add=False):
        import re

        with open(coords_input, 'rU') as infile:
            rawLinks = self.justLinks(infile.readlines());

        if algo == 'promer':
            shift = 4
        elif algo == 'nucmer':
            shift = 0

        contig2start_stop_list = {}

        link_file = coords_input.split('.')[0] + '.heat'
        alignment_list = []
        with open(link_file, 'w') as f:
            i = 1
            hit_list = []
            query_list = []

            for n, row in enumerate(rawLinks):

                l = re.split(r'\t+', row.rstrip('\n'))

                if l[9] not in contig2start_stop_list:
                    contig2start_stop_list[l[9+shift]] = {}
                    if contigs_add is not False:
                        contig2start_stop_list[l[9+shift]]["start"] = [int(l[0]) + contigs_add[l[9+shift]]]
                        contig2start_stop_list[l[9+shift]]["stop"] = [int(l[1]) + contigs_add[l[9+shift]]]
                    else:
                        contig2start_stop_list[l[9+shift]]["start"] = [int(l[0])]
                        contig2start_stop_list[l[9+shift]]["stop"] = [int(l[1])]
                else:
                    if contigs_add is not False:
                        contig2start_stop_list[l[9+shift]]["start"].append(int(l[0])+ contigs_add[l[9+shift]])
                        contig2start_stop_list[l[9+shift]]["stop"].append(int(l[1])+ contigs_add[l[9+shift]])
                    else:
                        contig2start_stop_list[l[9+shift]]["start"].append(int(l[0]))
                        contig2start_stop_list[l[9+shift]]["stop"].append(int(l[1]))


                # RhT_1 178 895 0
                contig = re.sub("\|", "", l[9+shift])
                if contigs_add:
                    start = int(l[0]) + contigs_add[contig]
                    end = int(l[1]) + contigs_add[contig]
                else:
                    start = int(l[0])
                    end = int(l[1])
                identity = l[6]
                alignment_list.append([contig, start, end, identity, identity])

                i += 1
                hit_list.append(re.sub("\|", "", l[9+shift]))
                query_list.append(re.sub("\|", "", l[10+shift]))
        return (alignment_list, hit_list, query_list, contig2start_stop_list)



    def megablast2heatmap(self, megablast_input):
        import re


        with open(megablast_input, 'rU') as infile:
            rawLinks = [i.rstrip().split('\t') for i in infile]

        contig2start_stop_list = {}

        heatmap_file = megablast_input.split('.')[0] + '.heat'
        with open(heatmap_file, 'w') as f:
            i = 1
            hit_list = []
            query_list = []
            for row in rawLinks:
                color = ''
                if row[0] not in contig2start_stop_list:
                    contig2start_stop_list[row[0]] = {}
                    contig2start_stop_list[row[0]]["start"] = [row[6]]
                    contig2start_stop_list[row[0]]["stop"] = [row[7]]
                else:
                    contig2start_stop_list[row[0]]["start"].append(row[6])
                    contig2start_stop_list[row[0]]["stop"].append(row[7])


                f.write(re.sub("\|", "", row[0]) + '\t' + row[6] + '\t' + row[7] + '\t' + row[2] + "\tz=%s\t" % row[2] +'\n')

                i += 1
                hit_list.append(re.sub("\|", "", row[0]))
                query_list.append(re.sub("\|", "", row[1]))
        return (hit_list, query_list, contig2start_stop_list)


    def get_gaps_from_nucmer_serie(self, fasta1, fasta_list):

        self.execute_promer(fasta1, fasta_list)

        all_gaps =  []

        for one_fasta in fasta_list:
            alignment_list, hit_list, query_list, contig2start_stop_list = self.nucmer_coords2data(one_fasta.split('.')[0] + '.coords')

            gaps = GapOverlapIdentification()
            one_nucmer_gaps = gaps.get_gap_data_from_start_stop_list(contig2start_stop_list, 0)

            all_gaps.append([one_fasta.split('.')[0], one_nucmer_gaps])
        return all_gaps

    def write_gap_files(self, gap_location_lists):
        for one_gap_list in gap_location_lists:
            for contig in one_gap_list[1]:
                with open( one_gap_list[0] + '.gaps', 'w') as f:
                    for one_gap in one_gap_list[1][contig]:
                        f.write('%s\t%s\t%s\n' % (contig, one_gap[0], one_gap[1]))

    def plot(self):
        '''
files <- dir('.', pattern=".gaps$")

y<-0

plot(c(0, 5000000), c(0, 250), type = "n", xlab = "", ylab = "")

for (i in 1:length(files)){
print(i)
table <- read.table(files[i], sep="\t", header=FALSE)
table[,4] <- table[,3]-table[,2]
table <- table[table[,4]>2000,]
apply(table, 1,function(x) rect(x[2],10+y,x[3],12+y, col="black"))
y <- y+4
}

all <- read.table("allgaps.tab", sep="\t", header=FALSE)
all[,4] <- all[,3]-all[,2]
all <- all[all[,4]>2000,]
r <- IRanges(all[,2], width =all[,4])
plot.new()
plotRanges(reduce(r))


apply(rr, 1,function(x) rect(x[2],10+y,x[3],12+y, col="black"))


plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                    col = "black", sep = 0.5, ...)
{
height <- 1
if (is(xlim, "Ranges"))
 xlim <- c(min(start(xlim)), max(end(xlim)))
bins <- disjointBins(IRanges(start(x), end(x) + 1))
plot.new()
plot.window(xlim, c(0, max(bins)*(height + sep)))
ybottom <- bins * (sep + height) - height
rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
title(main)
axis(1)
}
        '''



if __name__ == '__main__':
    ###Argument handling.
    arg_parser = argparse.ArgumentParser(description='');
    #arg_parser.add_argument("coords_input", help="Directory to show-coords tab-delimited input file.");
    arg_parser.add_argument("-r", "--fasta1", help="reference fasta")
    arg_parser.add_argument("-q", "--fasta2", help="query fasta", nargs='+')

    args = arg_parser.parse_args()

    comp = CompareWholeGenomes()
    comp.write_gap_files(comp.get_gaps_from_nucmer_serie(args.fasta1, list(args.fasta2)))

