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


class PHAST():
    def __init__(self, phast_database):

        self.phast_database = phast_database
        self.accession2classification = self.parse_phast_db()

    def parse_phast_db(self):
        import re
        search_patterns = {'transposase': re.compile(".*[tT]ransposase.*"),
                        'tail': re.compile(".*[tT]ail.*"),
                        'lysis': re.compile(".*[lL]ysis.*"),
                        'integrase': re.compile(".*[iI]ntegrase.*"),
                        'terminase': re.compile(".*[tT]erminase.*"),
                        'protease': re.compile(".*[pP]rotease.*"),
                        'fiber': re.compile(".*[fF]iber.*"),
                        'coat': re.compile(".*[cC]oat.*"),
                        'portal': re.compile(".*[pP]ortal.*")}

        accession2lassification = {}

        from Bio import SeqIO
        with open(self.phast_database) as f:
            records = SeqIO.parse(f, "fasta")
            for record in records:
                match = False
                for pattern in search_patterns:
                    if re.match(search_patterns[pattern], record.description) is not None:
                        #print pattern, record.description
                        accession2lassification[record.id] = pattern
                        match = True
                        break
                if match==False:
                    accession2lassification[record.id] = 'other'
        return accession2lassification

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]


class Genbank2circos():
    def __init__(self, genbank,
                 database,
                 samtools_depth=None,
                 gc=True,
                 accession2classification=False):


        import gbk2circos
        import shell_command


        self.fasta_aa = genbank.split(".")[0] + '.faa'
        self.fasta_nucl = genbank.split(".")[0] + '.fna'

        a,b,c = shell_command.shell_command("gbk2fna.py -i %s -o %s" % (genbank, self.fasta_nucl))
        a,b,c = shell_command.shell_command("gbk2faa.py -i %s -o %s -f" % (genbank, self.fasta_aa))

        print a,b,c

        self.circos_reference = gbk2circos.Circos_config("circos_contigs.txt", show_ideogram_labels="no", radius=0.5,show_tick_labels="yes", show_ticks="yes")

        self.get_karyotype_from_genbank(genbank)
        self.execute_blast(self.fasta_aa, database)

        class_colors = self.blast2barplot("blast_result.tab", bar_file="circos.bar", accession2classification=accession2classification)
        if not class_colors:
            histo_colors = 'blue'
            self.last_track = 0.7
        else:
            histo_colors = ''
            for class_col in class_colors:
                print class_col, class_colors[class_col]
                histo_colors+='%s,' % class_colors[class_col]
            histo_colors = histo_colors[0:-1]
            self.circos_reference.add_plot("best_hit.bar",type="histogram",r0="0.65r", r1="0.67r",color=histo_colors, fill_color=histo_colors, thickness=0)
            self.last_track = 0.65
        self.circos_reference.add_plot("circos.bar",type="histogram",r0="0.7r", r1="0.9r",color=histo_colors, fill_color=histo_colors, thickness=0)




        if samtools_depth is not None:
            for i, depth_file in enumerate(samtools_depth):
                self.samtools_depth2circos_data(depth_file, i)
                self.add_samtools_depth_track('circos_samtools_depth_%s.txt' % i)

        if gc:
            import GC
            from Bio import SeqIO

            fasta_records = list(SeqIO.parse(self.fasta_nucl, 'fasta'))

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

        self.config = self.circos_reference.get_file()



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



    def get_karyotype_from_genbank(self, genbank, out="circos_contigs.txt"):
        from Bio import SeqIO
        import re

        self.locus_tag2start_stop = {}

        genbank_data = [i for i in SeqIO.parse(open(genbank), "genbank")]

        #print "hit_list", hit_list
        with open(out, 'w') as f:
            # chr - Rhab Rhab 0 1879212 spectral-5-div-4
            i = 0
            contig_start = 0
            contig_end = 0
            for record in genbank_data:
                if i == 4:
                    i =0
                i+=1
                name = re.sub("\|", "", record.name)

                n4 = name
                if not 'n2' in locals():
                    n2 = name
                # spectral-5-div-%s
                f.write("chr - %s %s %s %s greys-3-seq-%s\n" % (name, name, 0, len(record.seq), i)) # , i

                for feature in record.features:
                    if feature.type == 'CDS':
                        self.locus_tag2start_stop[feature.qualifiers['locus_tag'][0]] = [record.name, int(feature.location.start), int(feature.location.end)]
        print self.locus_tag2start_stop

    def execute_blast(self,fasta, database):
        import shell_command

        cmd1 = "formatdb -i %s -p T" % database
        cmd2 = 'blastp -query %s -db %s -evalue 1e-5 -num_threads 8 -outfmt 6 -out blast_result.tab' % (fasta, database)
        #a, b, c = shell_command.shell_command(cmd1)
        a, b, c = shell_command.shell_command(cmd2)
        print a,b,c



    def blast2barplot(self, blast_output, bar_file="circos.bar", accession2classification=False):
        import re
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl

        if accession2classification:
            class2color = {}
            all_class = set(accession2classification.values())
            n_class = len(all_class)
            if n_class > 10:
                print all_class
                print n_class
                raise('Too mainy classes of match')
            else:
                for i, one_class in enumerate(all_class):
                    class2color[one_class] = 'paired-10-qual-%s' % (i+1)

        with open(blast_output, 'rU') as infile:
            infile = [i for i in infile]
            locus_list = [i.rstrip().split('\t')[0] for i in infile]
            hit_list = [i.rstrip().split('\t')[1] for i in infile]

        locus_tag2count = {}
        locus_tag2best_hit = {}
        for n, locus in enumerate(locus_list):
            if not accession2classification:
                if locus not in locus_tag2count:
                    locus_tag2count[locus] = 1
                else:

                    locus_tag2count[locus] += 1
            else:
                if locus not in locus_tag2count:
                    locus_tag2count[locus] = {}
                    locus_tag2best_hit[locus] = {}
                    for one_class in all_class:
                        locus_tag2count[locus][one_class] = 0
                        locus_tag2best_hit[locus][one_class] = 0
                    locus_tag2count[locus][accession2classification[hit_list[n]]] +=1
                    locus_tag2best_hit[locus][accession2classification[hit_list[n]]] +=1
                else:

                    locus_tag2count[locus][accession2classification[hit_list[n]]] +=1


        bar_file = bar_file
        best_hit_classif = open('best_hit.bar', 'w')
        with open(bar_file, 'w') as f:
            for locus in locus_tag2count:

                # RhT_1 178 895 0
                if not accession2classification:

                    f.write("%s\t%s\t%s\t%s\n" % (self.locus_tag2start_stop[locus][0],
                                                  self.locus_tag2start_stop[locus][1],
                                                  self.locus_tag2start_stop[locus][2],
                                                  locus_tag2count[locus]))
                else:
                    counts = [str(i) for i in locus_tag2count[locus].values()]
                    f.write("%s\t%s\t%s\t%s\n" % (self.locus_tag2start_stop[locus][0],
                                                  self.locus_tag2start_stop[locus][1],
                                                  self.locus_tag2start_stop[locus][2],
                                                  ','.join(counts)))
                    counts = [str(i) for i in locus_tag2best_hit[locus].values()]
                    best_hit_classif.write("%s\t%s\t%s\t%s\n" % (self.locus_tag2start_stop[locus][0],
                                                  self.locus_tag2start_stop[locus][1],
                                                  self.locus_tag2start_stop[locus][2],
                                                  ','.join(counts)))
        return class2color


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
    arg_parser.add_argument("-g", "--genbank", help="genbank")
    arg_parser.add_argument("-d", "--database", help="database")
    arg_parser.add_argument("-o", "--out_name", help="output name", default="circos_out")
    arg_parser.add_argument("-p", "--phast_database", action="store_true",help="PHAST database")
    args = arg_parser.parse_args()

    ##Run main
    if args.phast_database:
        p = PHAST(args.database)
        accession2classif = p.accession2classification

    else:
        accession2classif = False

    circosf = Genbank2circos(args.genbank, args.database, accession2classification=accession2classif)
    circosf.write_circos_files(circosf.config, circosf.brewer_conf)
    circosf.run_circos(out_prefix=args.out_name)



