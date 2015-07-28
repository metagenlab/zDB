#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# create html table from genome tabulated file
# headers: accession	size	gi	n proteins	n contigs	gc 	description
# add 4 columns with links of the form /assets/chlamdb/ffn/ for gbk/faa/ffm/fna
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2015
# ---------------------------------------------------------------------------


def tab2htm(input_file):
    with open(input_file, "r") as f:
        i = 0
        print '<table cellspacing="0" border="0">'

        for line in f:
            if i == 0:
                header = line.rstrip().split("\t")
                print '<tr>'
                print '    <th>%s</th>' % header[0]
                print '    <th>%s</th>' % header[1]
                print '    <th>%s</th>' % header[3]
                print '    <th>%s</th>' % header[4]
                print '    <th>%s</th>' % header[5]
                print '    <th>%s</th>' % header[6]
                print '    <th colspan="5">Download</th>'
                print '    </tr>'

                i+=1
                continue
            data = line.rstrip().split("\t")
            print '<tr>'
            print '    <td><a href="http://www.ncbi.nlm.nih.gov/nuccore/%s">%s<a></td>' % (data[2], data[0])
            print '    <td>%s</td>' % data[1]
            print '    <td>%s</td>' % data[3]
            print '    <td>%s</td>' % data[4]
            print '    <td>%s</td>' % data[5]
            print '    <td>%s</td>' % data[6]

            #print '<td colspan="6"><table width="800" border=0  class=table_genomes>'

            gbk = '<td style="vertical-align:top" nowrap>' \
            '<a  style="text-decoration:none;" href="/assets/chlamdb/gbk/%s.gbk">' \
            '<span style="background-color:#DCFFF0;color:black;padding:1px;">GBK</span></a>' \
            '</td>' % (data[0].split(".")[0])

            faa = '<td snowrap>' \
            '<a  style="text-decoration:none;" href="/assets/chlamdb/faa/%s.faa">' \
            '<span style="background-color:#DCFFF0;color:black;padding:1px;">FAA</span></a>' \
            '</td>' % (data[0].split(".")[0])

            fna = '<td nowrap>' \
            '<a  style="text-decoration:none;" href="/assets/chlamdb/fna/%s.fna">' \
            '<span style="background-color:#DCFFF0;color:black;padding:1px;">FNA</span></a>' \
            '</td>' % (data[0].split(".")[0])

            ffn = '<td nowrap>' \
            '<a  style="text-decoration:none;" href="/assets/chlamdb/ffn/%s.ffn">' \
            '<span style="background-color:#DCFFF0;color:black;padding:1px;">FFN</span></a>' \
            '</td>' % (data[0].split(".")[0])

            tab = '<td nowrap>' \
            '<a  style="text-decoration:none;" href="/assets/chlamdb/tab/%s.tab">' \
            '<span style="background-color:#DCFFF0;color:black;padding:1px;">TAB</span></a>' \
            '</td>' % (data[0].split(".")[0])

            #print '<tr>'
            print gbk
            print faa
            print fna
            print ffn
            print tab
            #print '</tr>'
            #print '</table></td>'
            #print '</tr>'

    print '</table>'

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_tab', type=str, help="input tabulated file")


    args = parser.parse_args()

    tab2htm(args.input_tab)


