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

def get_gaps_from_start_stop_lists(contig2start_stop_lists, contig2length=False, out_highlight="circos_gaps_highlight.txt", out_labels='circos_gaps_labels.txt'):
    import pandas as pd
    import numpy as np

    h = open(out_highlight, 'w')
    gap_data = {}
    for contig in contig2start_stop_lists:
        print contig
        #print contig2start_stop_lists[contig]["start"]
        data = pd.DataFrame({'start': contig2start_stop_lists[contig]["start"],
                'stop': contig2start_stop_lists[contig]["stop"] })
        data.start = data.start.astype(np.int64)
        data.stop = data.stop.astype(np.int64)
        #print 'add', data['start'][0]+1
        data_sort = data.sort(columns=["start"])
        #print data_sort
        index_start = 0
        comparison_index = 1
        stop = False

        while index_start < len(data_sort['start']):
            #print 'index', index_start
            #print 'compt', comparison_index
            if index_start+comparison_index >= len(data_sort['start']):
                break
            #print 'index', index_start
            #print 'comp index', comparison_index
            # cas si fragment suivant est compris a l'interieur du 1er
            # comparer avec le fragment suivant
            if int(data_sort['stop'][index_start+comparison_index])<= int(data_sort['stop'][index_start]):
                #print 'compris interne', index_start, comparison_index
                #print '1', data_sort['start'][index_start], data_sort['stop'][index_start]
                #print '2', data_sort['start'][index_start+comparison_index], data_sort['stop'][index_start+comparison_index]

                comparison_index+=1
                continue
            # si le start du 2eme se trouve dans le premier et est plus long
            # skipper le 1er fragment
            elif int(data_sort['start'][index_start+comparison_index])<=int(data_sort['stop'][index_start]) and  \
                int(data_sort['stop'][index_start+comparison_index]) > int(data_sort['stop'][index_start]):
                index_start+=comparison_index
                comparison_index=1
                continue
            elif int(data_sort['start'][index_start+comparison_index])-int(data_sort['stop'][index_start]) > 5000:
                print '############# Potential gap!!! ', contig, data_sort['stop'][index_start], "-" , data_sort['start'][index_start+comparison_index]
                if contig not in gap_data:
                    gap_data[contig] = [[data_sort['stop'][index_start], data_sort['start'][index_start+comparison_index]]]
                else:
                    gap_data[contig].append([data_sort['stop'][index_start], data_sort['start'][index_start+comparison_index]])


                #print index_start, comparison_index
                #print '1', data_sort['start'][index_start], data_sort['stop'][index_start]
                #print '2', data_sort['start'][index_start+comparison_index], data_sort['stop'][index_start+comparison_index]
                #print data_sort

                index_start+=comparison_index
                comparison_index=1
            else:
                #print 'add'
                index_start+=comparison_index
                comparison_index=1
    print "gap_data"
    print gap_data
    #import sys
    #sys.exit()

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
    #import sys
    #sys.exit()

'''
<plot>
	type		    = text
	 color              = black
	 file               = circos_gaps_labels.txt
	r0 = 1r
	r1 = 1.3r
	z = -100
	show_links     = yes
	link_dims      = 4p,4p,248p,4p,4p
	link_thickness = 2p
	link_color     = red

	label_size   = 24p
	label_font   = condensed

	padding  = 0p
	rpadding = 0p


 </plot>

'''


class Fasta2circos():
    def __init__(self, fasta1, fasta2,
                 filter_ref=True,
                 filter_query=True,
                 heatmap=False,
                 samtools_depth=False,
                 gaps=False,
                 blast=False,
                 gc=True,
                 highlight_list=[],
                 algo="nucmer"):

        print "fasta1", fasta1
        print "fasta2", len(fasta2), fasta2
        print "highlight list", highlight_list
        print "heatmap", heatmap
        if algo == "nucmer":
            self.execute_promer(fasta1, fasta2, algo="nucmer")
        elif algo == "megablast":
            self.execute_megablast(fasta1, fasta2)
        elif algo == "promer":
            self.execute_promer(fasta1, fasta2, algo="promer")
        if not heatmap:
            hit_list, query_list = self.get_link("%s.coords" % fasta2[0].split('.')[0])
        else:
            print 'wrinting link file'

            import gbk2circos
            self.circos_reference = gbk2circos.Circos_config("circos_contigs.txt", show_ideogram_labels="no", radius=0.5,show_tick_labels="yes", show_ticks="yes")
            #genome_list = [i.split('.')[0] + '.heat' for i in fasta2]

            #self.add_multiple_genome_tracks(genome_list, highlight_list)

            updated_list = []
            gap_hilight_list = []
            for i, one_fasta in enumerate(fasta2):
                if i%2 == 0:
                    col=200
                else:
                    col=250
                try:
                    if algo == "nucmer" or algo == 'promer':
                        hit_list, query_list, contig2start_stop_list = self.nucmer_coords2heatmap("%s.coords" % one_fasta.split('.')[0], col = col, algo=algo)
                    elif algo == "megablast":
                        hit_list, query_list, contig2start_stop_list = self.megablast2heatmap("blast_result_%s.tab" % one_fasta.split('.')[0], col = col)

                        raise IOError('unknown algo!')
                    updated_list.append(one_fasta)
                except:
                    continue
                if gaps:
                    print contig2start_stop_list
                    get_gaps_from_start_stop_lists(contig2start_stop_list, out_highlight="circos_gaps_highlight_%s.txt" % i, out_labels="circos_gaps_labels_%s.txt" % i)
                    if i%2 == 0:
                        color = "orrd-9-seq-9"
                    else:
                        color = "blues-9-seq-9"
                    gap_hilight_list.append([color, "circos_gaps_highlight_%s.txt" % i])


                    supp = '''show_links     = yes
                            link_dims      = 4p,4p,4p,4p,4p
                            link_thickness = 2p
                            link_color     = red

                            label_size   = 24p
                            label_font   = condensed

                            padding  = 0p
                            rpadding = 0p
                            '''
                    #self.circos_reference.add_plot("circos_gaps_labels_%s.txt" % i ,type="text", r0="1r", r1="1.3r",color="black",rules=supp)
            genome_list = [i.split('.')[0] + '.heat' for i in updated_list]
            self.add_multiple_genome_tracks(genome_list, highlight_list)
            if gaps:
                for i in gap_hilight_list:
                    self.circos_reference.add_highlight(i[1],
                                                        i[0],
                                                        r1="%sr" % (self.last_track - 0.01),
                                                        r0="%sr" % (self.last_track - 0.03))
                    self.last_track-=0.03



                    
            if blast:
                self.blast2circos_file(blast, fasta1)

                #self.circos_reference.add_highlight("circos_blast.txt",
                #                               'red',
                #                               r1="%sr" % (self.last_track - 0.01),
                #                               r0="%sr" % (self.last_track - 0.03))
                #self.last_track -= 0.03

                supp = '''
                        label_snuggle             = yes
                        snuggle_refine                 = yes
                        max_snuggle_distance            = 4r
                        show_links     = yes
                        link_dims      = 10p,8p,20p,4p,4p
                        link_thickness = 2p
                        link_color     = red

                        label_size   = 24p
                        label_font   = condensed

                        padding  = 0p
                        rpadding = 0p
                        '''
                self.circos_reference.add_plot("circos_blast_labels.txt", type="text", r0="1r", r1="1.3r",color="black", rules=supp)
            if gc:
                import GC
                from Bio import SeqIO

                fasta_records = list(SeqIO.parse(fasta1, 'fasta'))

                out_var_file = ('circos_GC_var.txt')
                out_skew_file = ('circos_GC_skew.txt')
                f = open(out_var_file, 'w')
                g = open(out_skew_file, 'w')

                out_var = ''
                out_skew = ''
                for record in fasta_records:
                    # this function handle scaffolds (split sequence when encountering NNNNN regions)
                    out_var += GC.circos_gc_var(record, 1000)

                    out_skew += GC.circos_gc_skew(record, 1000)
                #print out_skew
                f.write(out_var)
                g.write(out_skew)
                f.close()
                g.close()

                rule = """<rule>
                        condition          = var(value) < 0
                        fill_color         = lred
                        color = red
                        </rule>

                        <rule>
                        condition          = var(value) > 0
                        fill_color         = lblue
                        color = blue
                        </rule>
                """

                rule2 = """<rule>
                        condition          = var(value) < 0
                        fill_color         = lgreen
                        color = green
                        </rule>

                        <rule>
                        condition          = var(value) > 0
                        fill_color         = lblue
                        color = blue
                        </rule>
                """

                conditions = self.circos_reference.template_rules % (rule)
                self.circos_reference.add_plot('circos_GC_skew.txt', fill_color="green", r1="%sr" % (self.last_track -0.02), r0= "%sr" % (self.last_track -0.1), type="line", rules=conditions)
                conditions = self.circos_reference.template_rules % (rule2)
                self.circos_reference.add_plot('circos_GC_var.txt', fill_color="green", r1="%sr" % (self.last_track -0.12), r0= "%sr" % (self.last_track -0.2), type="line", rules=conditions)
                self.last_track = self.last_track -0.12

        if heatmap:
            c1, c2, c3, c4 = self.get_karyotype_from_fasta(fasta1, fasta2, hit_list, query_list, filter_ref, filter_query, both_fasta=False)
            #self.config = self.get_circos_config(c1, c2, c3, c4, link=False, heat=True)

            self.config = self.circos_reference.get_file()

            if samtools_depth is not None:
                for i, depth_file in enumerate(samtools_depth):
                    self.samtools_depth2circos_data(depth_file, i)
                    self.add_samtools_depth_track('circos_samtools_depth_%s.txt' % i)

        else:
            c1, c2, c3, c4 = self.get_karyotype_from_fasta(fasta1, fasta2, hit_list, query_list, filter_ref, filter_query, both_fasta=True)

            self.config = self.get_circos_config(c1, c2, c3, c4, heat=False)

        self.brewer_conf = """


        # Colors from www.ColorBrewer.org by Cynthia A. Brewer, Geography, Pennsylvania State University.
        # See BREWER for license. See www.colorbrewer.org for details.
        #
        # Color names are formatted as PALETTE-NUMCOLORS-TYPE-COLORCODE
        #
        # where PALETTE is the palette name, NUMCOLORS is the number of colors in the palette, TYPE is the palette type (div, seq, qual) and COLORCODE is the color index within the palette (another versison of the color is defined where COLORCODE is the color's letter, unique to a combination of PALETTE and TYPE)

        #
        # the value after the trailing comment is the palette's critical value. See http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html for details.

        # sequential palettes

        # 5 color palettes
        spectral-5-div = spectral-5-div-(\d+)
        spectral-5-div-rev = rev(spectral-5-div-(\d+))
        spectral-5-div-1 = 215,25,28 # 3
        spectral-5-div-c = 215,25,28 # 3
        spectral-5-div-2 = 253,174,97 # 3
        spectral-5-div-f = 253,174,97 # 3
        spectral-5-div-3 = 255,255,191 # 3
        spectral-5-div-h = 255,255,191 # 3
        spectral-5-div-4 = 171,221,164 # 3
        spectral-5-div-j = 171,221,164 # 3
        spectral-5-div-5 = 43,131,186 # 3
        spectral-5-div-m = 43,131,186 # 3
        """

    def blast2circos_file(self, blast, reference):

        import shell_command
        import best_blast_tree_display
        from Bio.Blast.Applications import NcbitblastnCommandline
        # todo catch IO errors, orther potential errors
        a,b,c = shell_command.shell_command('formatdb -i %s -p F' % (reference))
        #print a
        #print b
        print c
        tblastn_cline = NcbitblastnCommandline(query=blast,
                                         db=reference,
                                         evalue=0.001,
                                         outfmt=6,
                                         out="blast.tmp")

        stdout, stderr = tblastn_cline()

        #a,b,c = shell_command.shell_command('tblastn -query %s -db %s -evalue 1e-5 -max_target_seqs 1 -outfmt 6 > blast.tmp' % (blast, reference))
        #a,b,c = shell_command.shell_command('tblastn -query %s -db %s -evalue 1e-5 -max_target_seqs 1 -outfmt 6' % (blast, reference))
        print '############## BLAST ###################'
        #print a
        #print b
        #print c

        blast2data, queries = best_blast_tree_display.remove_blast_redundancy(["blast.tmp"])
        print blast2data


        o = open('circos_blast.txt', "w")
        l = open('circos_blast_labels.txt', "w")

        #with open(blast, 'r') as b:
        '''
        for line in a.split('\n'):

            data = line.rstrip().split('\t')
            #print data
            try:
                if float(data[2])>80:
                    location = sorted([data[8], data[9]])
                    o.write("%s\t%s\t%s\n" % (data[1], location[0], location[1]))
                    l.write("%s\t%s\t%s\t%s\n" % (data[1],  location[0], location[1], data[0]))
            except IndexError:
                continue
        '''
        for contig in blast2data:
            for gene in blast2data[contig]:
                if float(blast2data[contig][gene][0])>40:
                    location = sorted(blast2data[contig][gene][1:3])
                    o.write("%s\t%s\t%s\n" % (contig, location[0], location[1]))
                    l.write("%s\t%s\t%s\t%s\n" % (contig,  location[0], location[1], gene))

        o.close()

    def add_multiple_genome_tracks(self, track_file_list, highlight_list=[]):
        print 'track file list', track_file_list
        import os
        r1 = 0.95
        r0 = 0.93
        n = 0
        for i, orthofile in enumerate(track_file_list):
            n+=1
            if orthofile not in highlight_list:
                #print orthofile
                if n%2==0: # orrd-9-seq # blues
                    self.circos_reference.add_plot(orthofile, type="heatmap", r1="%sr" % r1, r0= "%sr" % r0, color="ylorrd-9-seq", fill_color="", thickness = "2p", z = 1, rules ="", backgrounds="",url="")
                else:
                    self.circos_reference.add_plot(orthofile, type="heatmap", r1="%sr" % r1, r0= "%sr" % r0, color="gnbu-9-seq", fill_color="", thickness = "2p", z = 1, rules ="", backgrounds="",url="")

                #circos.add_highlight(orthofile, fill_color="ortho3", r1="%sr" % r1, r0= "%sr" % r0, href=href)
                r1 = r1-0.023 # 046
                r0 = r0-0.023 # 046
            else:
                n-=1
                self.circos_reference.add_highlight(orthofile, r1="0.99r", r0= "0.97r",fill_color="piyg-9-div-8")
        self.last_track = r0
        self.config = self.circos_reference.get_file()

        #t = open('circos.config', "w")


    def add_samtools_depth_track(self, samtools_file):

        rules = """
        <rules>
        <rule>
        condition          = var(value) > 500
        show               = no
        </rule>

        <rule>
        condition          = var(value) < 100
        color              = red
        </rule>
        </rules>

        """

        self.circos_reference.add_plot(samtools_file,
                        type="histogram",
                        r1="%sr" % (self.last_track - 0.08),
                        r0= "%sr" % (self.last_track - 0.15),
                        color="black",
                        fill_color="grey_a5",
                        thickness = "1p",
                        z = 1,
                        rules =rules,
                        backgrounds="",
                        url="")
        self.last_track -= 0.15

        self.config = self.circos_reference.get_file()


    def run_circos(self, config_file="circos.config", out_prefix="circos"):
        import shell_command
        cmd = 'circos -outputfile %s.svg -conf %s' % (out_prefix, config_file)
        a,b,c = shell_command.shell_command(cmd)
        print "out", a, "err", b, "code", c
        if c == 255:
            raise CircosException("Circos problem, check files... quitting")


    def clean_tmp_files(self):
        import os
        os.remove("colors.my")
        os.remove("out.delta")
        #os.remove("circos.link")
        #os.remove("circos.html")
        os.remove("circos_contigs.txt")
        #os.remove("circos.config")
        os.remove("brewer.all.conf")
        d = os.getcwd()
        purge(d, ".*.heat")
        purge(d, "circos_gaps_labels.*")
        purge(d, "circos_gaps_highlight.*")
        purge(d, ".*.coords")
        purge(d, ".*a.n.*")

    def write_circos_files(self, config, color_config):

        with open("circos.config", 'w') as f:
            f.write(config)
        with open("brewer.all.conf", "w") as f:
            f.write(color_config)

    def id_generator(self, size=6, chars=string.ascii_lowercase ): # + string.digits
        return ''.join(random.choice(chars) for _ in range(size))


    def get_circos_config(self, last_seq_id, first_seq_id, mid1, mid2, link=True, heat=False):
        print 'config heatmap', heat

        link_code = '''

        <links>

        <link>
        ribbon = yes
        file          = circos.link
        color         = orange
        radius        = 0.95r
        bezier_radius = 0.1r
        thickness     = 3

        #<rules>
        #<rule>
        # you can also test whether only one end is
        # reversed using var(inv)
        #condition  = var(rev1) && ! var(rev2)
        #color      = blue
        #</rule>
        #<rule>
        #condition  = var(rev2) && ! var(rev1)
        #color      = blue
        #</rule>
        #<rule>
        #condition  = var(rev1) && var(rev2)
        #color      = orange
        #</rule>

        #</rules>

        </link>

        </links>


        '''

        plot_templat = '''

        <plots>

        %s

        </plots>



        '''



        heatmap = '''
        <plot>
        type		    = heatmap
         r0                 = 0.9r
         r1                 = 0.97r
         color              = blues-6-seq-3, orrd-9-seq
         fill_color         =
         thickness          = 2p
         file               = circos.heat
         z                  = 1

         </plot>
        '''


        chr_spacing = '''

           <pairwise %s %s>
         spacing = 2u
        </pairwise>

        '''


        circos_config = '''

         karyotype = circos_contigs.txt
         chromosomes_units           = 10000
         chromosomes_display_default = yes
        <ideogram>

         <spacing>
         default            = %su

         %s

         </spacing>



         # thickness and color of ideograms
         thickness          = 15p
         stroke_thickness   = 1
         stroke_color       = black

         # the default chromosome color is set here and any value
         # defined in the karyotype file overrides it
         fill               = yes
         fill_color         = black

         # fractional radius position of chromosome ideogram within image
         radius             = 0.65r
         show_label         = no
         label_font         = default
         label_radius       = dims(ideogram,radius) + 0.175r
         label_size         = 30
         label_parallel     = no

         # show_bands determines whether the outline of cytogenetic bands
         # will be seen
         show_bands         = yes
         band_stroke_thickness = 1
          # in order to fill the bands with the color defined in the karyotype
         # file you must set fill_bands
         fill_bands         = yes
         band_transparency  = 1

         </ideogram>
        show_ticks         = yes
         show_tick_labels   = yes

         <ticks>

         tick_label_font    = condensed
         radius             = dims(ideogram,radius_outer)
         label_offset       = 8p
         label_size         = 8p
         color              = black
         thickness          = 4p

         #<tick>
         #spacing           = 100u
         #size              = 16p
         #label_multiplier  = 1e-3
         #show_label        = yes
         #label_size        = 35p
         #format            = %%d kb
         #thickness         = 5p
         #</tick>


         <tick>
         skip_first_label = yes
         multiplier   = 1/1u
         spacing           = 1u
         size              = 15p
         show_label        = no
         label_size        = 25p
         </tick>

         <tick>
         skip_first_label = yes
         multiplier   = 10/1u
         spacing           = 10u
         size              = 10p
         show_label        = yes
         label_size        = 25p
         format            = %%.1d kb
         </tick>

        <tick>
        multiplier   = 1
        position = end
        size              = 16p
        show_label        = yes
        label_size        = 25p
        format    = %%d bp
        </tick>

        <tick>
        position = 0u
        size              = 16p
        label_size        = 25p
        label = a
        </tick>

         </ticks>


        %s



        <colors>
          <<include colors.my>>
         <<include brewer.all.conf>>
         </colors>
         <image>
         image_map_use      = yes
         image_map_overlay  = yes
         image_map_overlay_stroke_color     = red
         <<include etc/image.conf>>
         </image>
         <<include etc/colors_fonts_patterns.conf>>
         # system and debug settings
         <<include etc/housekeeping.conf>>
         anti_aliasing*     = no

         '''

        if link:

            if last_seq_id == mid1 and mid2 == first_seq_id:
                ch = chr_spacing % (last_seq_id, first_seq_id)
                return circos_config % (0.5, ch, link_code)
            else:
                ch1 = chr_spacing % (last_seq_id, first_seq_id)
                ch2 = chr_spacing % (mid1, mid2)
                return circos_config % (0.5, ch1 + ch2, link_code)

        if heat:
            print '############ ok ###########'
            return circos_config % (0, "", plot_template % heatmap)












    def get_karyotype_from_fasta(self, fasta1, fasta2, hit_list, query_list, filter_ref=True, filter_query=True, out="circos_contigs.txt", both_fasta = True):
        print 'both fasta', both_fasta
        from Bio import SeqIO
        import re
        print 'hit list', hit_list
        fasta_data1 = [i for i in SeqIO.parse(open(fasta1), "fasta")]
        fasta_data2 = [i for i in SeqIO.parse(open(fasta2[0]), "fasta")]
        #print "hit_list", hit_list
        with open(out, 'w') as f:
            # chr - Rhab Rhab 0 1879212 spectral-5-div-4
            i = 0
            contig_start = 0
            contig_end = 0
            for record in fasta_data1:
                if i == 4:
                    i =0
                i+=1
                name = re.sub("\|", "", record.name)
                if filter_ref:
                    if name in hit_list:
                        n4 = name
                        if not 'n2' in locals():
                            n2 = name
                        # spectral-5-div-%s
                        f.write("chr - %s %s %s %s greys-3-seq-%s\n" % (name, name, 0, len(record), i)) # i,
                else:
                    n4 = name
                    if not 'n2' in locals():
                        n2 = name
                    # spectral-5-div-%s
                    f.write("chr - %s %s %s %s greys-3-seq-%s\n" % (name, name, 0, len(record), i)) # , i

            if both_fasta:
                for record in fasta_data2:
                    name = re.sub("\|", "", record.name)
                    if filter_query:
                            if name in query_list:
                                f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (name, name, 0, len(record), i))
                    else:
                        f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (name, name, 0, len(record), i))
        n1 = re.sub("\|", "", fasta_data2[-1].name)
        n3 = re.sub("\|", "", fasta_data2[0].name)
        if not 'n2' in locals():
            n2 = n1
        return (n1, n2, n3, n4)

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


    def get_link(self, coords_input, link_file="circos.link"):
        import re
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl


        with open(coords_input, 'rU') as infile:
            rawLinks = self.justLinks(infile.readlines());

        all_id = [float(re.split(r'\t+', i.rstrip('\n'))[6]) for i in rawLinks]
        all_start = []
        all_ends = []
        #print len(all_id), all_id
        #print sorted(all_id)

        id_min = min(all_id)
        id_max = max(all_id)
        #print 'mimax', id_min, id_max
        norm = mpl.colors.Normalize(vmin=id_min, vmax=id_max)
        cmap = cm.Blues
        cmap_blue = cm.OrRd
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
        c = open('colors.my', 'w')

        with open(link_file, 'w') as f:
            #c.write()
            i = 1
            hit_list = []
            query_list = []
            for row in rawLinks:
                color = ''
                l = re.split(r'\t+', row.rstrip('\n'))

                # print l[14]
                #f.write(l[13] + '\t' + l[0] + '\t' + l[1] + '\t' + re.sub("\|", "", l[14]) + '\t' + l[2] + '\t' + l[3] + '\n')
                if int(l[0]) < int(l[1]) and int(l[2]) > int(l[3]):
                    color = m.to_rgba(float(l[6]))

                elif int(l[0]) > int(l[1]) and int(l[2]) < int(l[3]):
                    color = m.to_rgba(float(l[6]))

                else:
                    color = m2.to_rgba(float(l[6]))
                #print color, float(l[6])
                color_id = self.id_generator()
                #print color_id, 'id:',l[6]
                c.write('%s = %s,%s,%s,%s\n' % (color_id,
                                                int(round(color[0]*250,0)),
                                                int(round(color[1]*250,0)),
                                                int(round(color[2]*250,0)),
                                                0.5))

                f.write(re.sub("\|", "", l[9]) + '\t' + l[0] + '\t' + l[1] + '\t' + re.sub("\|", "", l[10]) + '\t' + l[2] + '\t' + l[3] + '\tcolor=%s' % color_id +'\n')
                # sys.stdout.write
                i += 1
                hit_list.append(re.sub("\|", "", l[9]))
                query_list.append(re.sub("\|", "", l[10]))
        return (hit_list, query_list)



    def nucmer_coords2heatmap(self, coords_input, link_file="circos.heat", col=250, algo="nucmer"):
        import re
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl


        with open(coords_input, 'rU') as infile:
            rawLinks = self.justLinks(infile.readlines());

        all_id = [float(re.split(r'\t+', i.rstrip('\n'))[6]) for i in rawLinks]

        #print len(all_id), all_id

        #print sorted(all_id)
        try:
            id_min = min(all_id)
        except:
            return None
        id_max = max(all_id)
        #print 'mimax', id_min, id_max
        norm = mpl.colors.Normalize(vmin=id_min, vmax=id_max)
        cmap = cm.Blues
        cmap_blue = cm.OrRd
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
        c = open('colors.my', 'w')


        if algo == 'promer':
            shift = 4
        elif algo == 'nucmer':
            shift = 0

        contig2start_stop_list = {}

        link_file = coords_input.split('.')[0] + '.heat'
        with open(link_file, 'w') as f:
            i = 1
            hit_list = []
            query_list = []

            for n, row in enumerate(rawLinks):
                color = ''
                l = re.split(r'\t+', row.rstrip('\n'))
                if n == 0:
                    f.write("%s\t%s\t%s\t0\tz=0\n" % (re.sub("\|", "", l[9+shift]), l[0], int(l[0])+1))
                    f.write("%s\t%s\t%s\t100\tz=0\n" % (re.sub("\|", "", l[9+shift]), l[0], int(l[0])+1))
                if l[9] not in contig2start_stop_list:
                    contig2start_stop_list[l[9+shift]] = {}
                    contig2start_stop_list[l[9+shift]]["start"] = [l[0]]
                    contig2start_stop_list[l[9+shift]]["stop"] = [l[1]]
                else:
                    contig2start_stop_list[l[9+shift]]["start"].append(l[0])
                    contig2start_stop_list[l[9+shift]]["stop"].append(l[1])
                # print l[14]
                #f.write(l[13] + '\t' + l[0] + '\t' + l[1] + '\t' + re.sub("\|", "", l[14]) + '\t' + l[2] + '\t' + l[3] + '\n')
                if int(l[0]) < int(l[1]) and int(l[2]) > int(l[3]):
                    color = m.to_rgba(float(l[6]))

                elif int(l[0]) > int(l[1]) and int(l[2]) < int(l[3]):
                    color = m.to_rgba(float(l[6]))

                else:
                    color = m2.to_rgba(float(l[6]))
                #print color, float(l[6])
                color_id = self.id_generator()
                #print color_id, 'id:',l[6]

                c.write('%s = %s,%s,%s,%s\n' % (color_id,
                                                int(round(color[0]*col,0)),
                                                int(round(color[1]*col,0)),
                                                int(round(color[2]*col,0)),
                                                0.5))

                # RhT_1 178 895 0

                f.write(re.sub("\|", "", l[9+shift]) + '\t' + l[0] + '\t' + l[1] + '\t' + l[6] + "\tz=%s\t" % l[6] +'\n')
                # sys.stdout.write
                i += 1
                hit_list.append(re.sub("\|", "", l[9+shift]))
                query_list.append(re.sub("\|", "", l[10+shift]))
        return (hit_list, query_list, contig2start_stop_list)


    def megablast2heatmap(self, megablast_input, link_file="circos.heat", col=250):
        import re
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl


        with open(megablast_input, 'rU') as infile:
            rawLinks = [i.rstrip().split('\t') for i in infile]

        all_id = [float(i[2]) for i in rawLinks]

        #print len(all_id), all_id

        #print sorted(all_id)
        try:
            id_min = min(all_id)
        except:
            return None
        id_max = max(all_id)
        #print 'mimax', id_min, id_max
        norm = mpl.colors.Normalize(vmin=id_min, vmax=id_max)
        cmap = cm.Blues
        cmap_blue = cm.OrRd
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
        c = open('colors.my', 'w')


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
                # print l[14]
                #f.write(l[13] + '\t' + l[0] + '\t' + l[1] + '\t' + re.sub("\|", "", l[14]) + '\t' + l[2] + '\t' + l[3] + '\n')
                color = m2.to_rgba(float(row[2]))
                #print color, float(l[6])
                color_id = self.id_generator()
                #print color_id, 'id:',l[6]

                c.write('%s = %s,%s,%s,%s\n' % (color_id,
                                                int(round(color[0]*col,0)),
                                                int(round(color[1]*col,0)),
                                                int(round(color[2]*col,0)),
                                                0.5))

                # RhT_1 178 895 0

                f.write(re.sub("\|", "", row[0]) + '\t' + row[6] + '\t' + row[7] + '\t' + row[2] + "\tz=%s\t" % row[2] +'\n')
                # sys.stdout.write
                i += 1
                hit_list.append(re.sub("\|", "", row[0]))
                query_list.append(re.sub("\|", "", row[1]))
        return (hit_list, query_list, contig2start_stop_list)


    def samtools_depth2circos_data(self, samtools_depth_file, i):
        import numpy
        contig2coverage = {}
        with open(samtools_depth_file, 'r') as f:
            for line in f:
                data = line.rstrip().split('\t')
                if data[0] not in contig2coverage:
                    contig2coverage[data[0]] = []
                    contig2coverage[data[0]].append(int(data[2]))
                else:
                    contig2coverage[data[0]].append(int(data[2]))
        with open('circos_samtools_depth_%s.txt' % i, 'w') as g:
            for contig in contig2coverage:
                # split list by chunks of 1000
                mychunks = chunks(contig2coverage[contig], 1000)
                for i, cov_list in enumerate(mychunks):
                    #print cov_list
                    median_depth = numpy.median(cov_list)
                    g.write("%s\t%s\t%s\t%s\n" % (contig, i*1000, (i*1000)+999, median_depth))



if __name__ == '__main__':
    ###Argument handling.
    arg_parser = argparse.ArgumentParser(description='');
    #arg_parser.add_argument("coords_input", help="Directory to show-coords tab-delimited input file.");
    arg_parser.add_argument("-r", "--fasta1", help="reference fasta")
    arg_parser.add_argument("-q", "--fasta2", help="query fasta", nargs='+')
    arg_parser.add_argument("-fr", "--filterr", action="store_false", help="do not remove reference sequences without any similarity from the plot (default False)")
    arg_parser.add_argument("-fq", "--filterq", action="store_false", help="do not remove query sequences without any similarity from the plot (default False)")
    arg_parser.add_argument("-l", "--link", help="link circos and not heatmap circos", action="store_true")
    arg_parser.add_argument("-s", "--samtools_depth", help="samtools depth file", nargs="+")
    arg_parser.add_argument("-o", "--output_name", help="output circos pefix", default="nucmer2circos")
    arg_parser.add_argument("-g", "--gaps", help="highlight gaps", action="store_true")
    arg_parser.add_argument("-b", "--blast", help="highlight blast hits (-outfmt 6)")
    arg_parser.add_argument("-n", "--highlight", help="highlight instead of heatmap corresponding list of records", nargs="+")
    arg_parser.add_argument("-a", "--algo", help="algorythm to use to compare the genome (megablast, nucmer or promer)", default="nucmer")
    args = arg_parser.parse_args()

    if args.highlight is None:
        args.highlight=[]
    ###Variable Definitions

    ##Run main
    circosf = Fasta2circos(args.fasta1,
                           args.fasta2,
                           args.filterr,
                           args.filterq,
                           heatmap=args.link,
                           samtools_depth=args.samtools_depth,
                           gaps=args.gaps,
                           blast=args.blast,
                           highlight_list=args.highlight,
                           algo=args.algo)

    circosf.write_circos_files(circosf.config, circosf.brewer_conf)
    circosf.run_circos(out_prefix=args.output_name)
    #circosf.clean_tmp_files()
