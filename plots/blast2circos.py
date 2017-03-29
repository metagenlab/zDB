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
import circos_utils

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
                 tabulated_blast_file = False,
                 samtools_depth=None,
                 gc=True,
                 accession2classification=False,
                 identity_cutoff=50,
                 query_coverage_cutoff=0.6,
                 evalue_cutoff=0.00005,
                 execute_blast=True,
                 blast_tab_file=False):


        import gbk2circos
        import shell_command
        from Bio import SeqIO



        self.reference_records = [i for i in SeqIO.parse(open(genbank), 'genbank')]
        self.contig_list = [record.name for record in self.reference_records]
        self.contig2gc_percent = circos_utils.contig2gc(self.reference_records)

        if samtools_depth:
            self.contig2median_depth = circos_utils.get_median_contig_coverage(self.contig_list, samtools_depth)
        else:
            self.contig2median_depth = False

        self.records2locus_tag2start_stop = circos_utils.records2locus_tag2start_stop(self.reference_records)
        # coordinates to get concatenated chromosome coord
        self.contigs_add = circos_utils.get_contigs_coords(self.reference_records)


        if not blast_tab_file:
            self.blast_tab_file = "blast_result.tab"
        else:
            self.blast_tab_file = blast_tab_file
        self.fasta_aa = genbank.split(".")[0] + '.faa'
        self.fasta_nucl = genbank.split(".")[0] + '.fna'
        print 'extractint aa and fna sequences...'
        a,b,c = shell_command.shell_command("gbk2fna.py -i %s -o %s" % (genbank, self.fasta_nucl))
        a,b,c = shell_command.shell_command("gbk2faa.py -i %s -f -o %s" % (genbank, self.fasta_aa))

        self.fasta_aa_records = [i for i in SeqIO.parse(open(self.fasta_aa), 'fasta')]
        self.fasta_nucl_records = [i for i in SeqIO.parse(open(self.fasta_nucl), 'fasta')]

        #else:
        #    self.blast_tab_file = tabulated_blast_file

        if execute_blast:
            print 'executing BLASTP...'
            self.execute_blast(self.fasta_aa, database)

        self.circos_reference = gbk2circos.Circos_config("circos_contigs.txt",
                                                         show_ideogram_labels="yes",
                                                         radius=0.7,
                                                         show_tick_labels="yes",
                                                         show_ticks="yes")


        self.first_contig, self.last_contig = circos_utils.get_karyotype_from_gbk_or_fasta_single_ref(self.reference_records,
                                                                                                      out="circos_contigs.txt")

        class_colors = self.blast2barplot(self.fasta_aa_records,
                                          database,
                                          self.blast_tab_file,
                                          bar_file="circos.bar",
                                          accession2classification=accession2classification,
                                          identity_cutoff=identity_cutoff,
                                          query_coverage_cutoff=query_coverage_cutoff,
                                          evalue_cutoff=evalue_cutoff)

        countour_col = 'black'

        if not class_colors:
            histo_colors = 'orange'

            self.last_track = 0.7
        else:
            histo_colors = ''

            for class_col in class_colors:
                print class_col, class_colors[class_col]
                histo_colors+='%s,' % class_colors[class_col]
            histo_colors = histo_colors[0:-1]
            self.circos_reference.add_plot("best_hit.bar",
                                           type="histogram",
                                           r0="0.65r", r1="0.67r",
                                           color=histo_colors,
                                           fill_color=histo_colors,
                                           thickness=0,
                                           min=0,
                                           max=50)
            self.last_track = 0.65

        add = '''
                <axes>
                <axis>
                spacing   = 0.1r
                color     = lgrey
                thickness = 2
                </axis>
                </axes>
                '''
        self.circos_reference.add_plot("circos.bar",
                                       type="histogram",
                                       r0="0.7r",
                                       r1="0.9r",
                                       color=countour_col,
                                       fill_color=histo_colors,
                                       thickness=1,
                                       min="0",
                                       max="50",
                                       rules=add)

        if samtools_depth:
            #for i, depth_file in enumerate(samtools_depth):
            all_contigs_median = circos_utils.samtools_depth2circos_data(samtools_depth,
                                                    self.contigs_add,
                                                    1)
            self.add_samtools_depth_track('circos_samtools_depth_1.txt',
                                          lower_cutoff=int(all_contigs_median)/2,
                                          top_cutoff=int(all_contigs_median)*2)

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
                out_var += GC.circos_gc_var(record,
                                            1000,
                                            shift=self.contigs_add[record.name][0])

                out_skew += GC.circos_gc_skew(record,
                                              1000,
                                              shift=self.contigs_add[record.name][0])
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



    def add_samtools_depth_track(self, samtools_file, lower_cutoff=50, top_cutoff=5000):

        rules = """
        <rules>
        <rule>
        condition          = var(value) > %s
        color              = green
        fill_color         = lgreen
        </rule>

        <rule>
        condition          = var(value) < %s
        color              = red
        fill_color         = lred
        </rule>
        </rules>

        """ % (top_cutoff,
               lower_cutoff)

        self.circos_reference.add_plot(samtools_file,
                        type="histogram",
                        r1="%sr" % (self.last_track - 0.08),
                        r0= "%sr" % (self.last_track - 0.20),
                        color="black",
                        fill_color="grey_a5",
                        thickness = "1p",
                        z = 1,
                        rules =rules,
                        backgrounds="",
                        url="")
        self.last_track -= 0.20

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



       #print self.locus_tag2start_stop

    def execute_blast(self,fasta, database):
        import shell_command

        cmd1 = "formatdb -i %s -p T" % database
        cmd2 = 'blastp -query %s -db %s -evalue 1e-5 -num_threads 8 -outfmt 6 -out blast_result.tab' % (fasta, database)
        #a, b, c = shell_command.shell_command(cmd1)
        a, b, c = shell_command.shell_command(cmd2)
        print a,b,c

    def read_blast_tab_file(self, query_records, input_database, blast_file, identity_cutoff, evalue_cutoff, query_cov_cutoff):
        from Bio import SeqIO
        # get the length of each query sequence
        # calculate coverage and filter hits
        locus2length = {}
        for record in query_records:
            locus2length[record.id] = len(record.seq)


        with open(input_database, 'r') as f:
            accession2description = {}
            records = SeqIO.parse(f,'fasta')
            for record in records:
                accession2description[record.id] = record.description

        with open(blast_file, 'rU') as infile:
            infile = [i.rstrip().split('\t') for i in infile]
            # KpGe_00001	gi|410609158|ref|YP_006952151.1|	82.88	146	25	0	1	146	1	146	4e-85	 250
            locus_list = []
            hit_list = []
            keep = []
            for line in infile:
                query_name = line[0]
                hit_name = line[1]
                query_align_length = (int(line[7])-int(line[6]))+1
                query_cov = query_align_length/float(locus2length[query_name])
                e_value = float(line[10])
                identity = float(line[2])
                if e_value <= evalue_cutoff and identity>=identity_cutoff and query_cov >= query_cov_cutoff:
                    keep.append(line)
                    hit_list.append(hit_name)
                    locus_list.append(query_name)

        with open('blast_best_hit.txt', 'w') as o:
            for n,line in enumerate(keep):
                try:
                    acc = accession2description[line[1]].split('| ')[1]
                except:
                    try:
                        acc = accession2description[line[1]].split('|')[1]
                    except:
                        acc='-'
                if n==0:
                    locus = line[0]
                    o.write('\t'.join(line) + '\t%s\n' % acc)
                else:
                    if locus == line[0]:
                        continue
                    else:
                        locus=line[0]

                        o.write('\t'.join(line) + '\t%s\n' % acc)

        return locus_list, hit_list

    def blast2barplot(self, input_records,
                      input_database,
                      blast_output,
                      bar_file="circos.bar",
                      accession2classification=False,
                      identity_cutoff=50,
                      query_coverage_cutoff=0.6,
                      evalue_cutoff=0.00005):
        import re
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl

        # group blast into categories
        # categories will be coloured differently on the circos figure
        if accession2classification:
            class2color = {}
            all_class = set(accession2classification.values())
            n_class = len(all_class)
            if n_class > 10:
                #print all_classblast2barplot
                #print n_class
                raise('Too mainy classes of match')
            else:
                for i, one_class in enumerate(all_class):
                    class2color[one_class] = 'paired-10-qual-%s' % (i+1)
        else:
            class2color = False

        # parse blast
        '''
        with open(blast_output, 'rU') as infile:
            infile = [i for i in infile]
            locus_list = [i.rstrip().split('\t')[0] for i in infile]
            hit_list = [i.rstrip().split('\t')[1] for i in infile]
        '''
        print 'identity', identity_cutoff
        print 'cov', query_coverage_cutoff
        print 'evalue', evalue_cutoff
        locus_list, hit_list = self.read_blast_tab_file(query_records=input_records,
                                                        blast_file=blast_output,
                                                        input_database=input_database,
                                                        identity_cutoff=identity_cutoff, # 80
                                                        evalue_cutoff=evalue_cutoff,
                                                        query_cov_cutoff=query_coverage_cutoff)

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
        contig_counts = open('coutig_couts.tab', 'w')
        if not self.contig2median_depth:
            contig_counts.write("contig\tN_blast\tN_no_blast\tratio\tGC\n")
        else:
            contig_counts.write("contig\tN_blast\tN_no_blast\tratio\tGC\tmedian_depth\tassembly_depth\tratio_depth\n")
        with open(bar_file, 'w') as f:
                for record in self.reference_records:
                    contig_hit_count = 0
                    contig_cds_count = 0
                    for feature in record.features:
                        if feature.type == 'CDS':
                            start = int(feature.location.start)
                            stop = int(feature.location.end)
                            locus = feature.qualifiers['locus_tag'][0]
                            try:
                                count = locus_tag2count[locus]
                                # set to 50 of to much hits
                                if count > 50:
                                    count = 50
                            except KeyError:
                                count = 0
                            if count > 0:
                                contig_hit_count+=1
                            if not 'pseudo' in feature.qualifiers:
                                contig_cds_count+=1
                            # RhT_1 178 895 0
                            if not accession2classification:
                                print start
                                print stop
                                print self.contigs_add[record.name]
                                f.write("%s\t%s\t%s\t%s\tid=%s\n" % (record.name,
                                                              str(int(start) + self.contigs_add[record.name][0]),
                                                              str(int(stop) + self.contigs_add[record.name][0]),
                                                              count,
                                                              locus))
                            else:
                                counts = [str(i) for i in locus_tag2count[locus].values()]
                                f.write("%s\t%s\t%s\t%s\n" % (record.name,
                                                              str(self.records2locus_tag2start_stop[locus][0] + self.contigs_add[record.name][0]),
                                                              str(self.records2locus_tag2start_stop[locus][1] + self.contigs_add[record.name][0]),
                                                              ','.join(counts)))
                                counts = [str(i) for i in locus_tag2best_hit[locus].values()]
                                best_hit_classif.write("%s\t%s\t%s\t%s\n" % (record.name,
                                                              str(self.records2locus_tag2start_stop[locus][0] + self.contigs_add[record.name][0]),
                                                              str(self.records2locus_tag2start_stop[locus][1] + self.contigs_add[record.name][0]),
                                                              ','.join(counts)))

                    print 'hits', contig_hit_count
                    print 'cds', contig_cds_count
                    if contig_cds_count == 0 and contig_hit_count == 0:
                        ratio = 0
                    else:
                        ratio = round(contig_hit_count/float(contig_cds_count), 2)
                    if not self.contig2median_depth:
                        contig_counts.write("%s\t%s\t%s\t%s\t%s\n" % (record.name,
                                                                  contig_cds_count,
                                                                  contig_hit_count,
                                                                  round(self.contig2gc_percent[record.name],2),
                                                                  ratio))
                    else:
                        ratio_depth = round(self.contig2median_depth[record.name]/float(self.contig2median_depth["assembly"]),2)
                        contig_counts.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (record.name,
                                                                                  contig_cds_count,
                                                                                  contig_hit_count,
                                                                                  ratio,
                                                                                  round(self.contig2gc_percent[record.name],2),
                                                                                  self.contig2median_depth[record.name],
                                                                                  self.contig2median_depth["assembly"],
                                                                                  ratio_depth))
        contig_counts.close()
        return class2color


if __name__ == '__main__':
    ###Argument handling.
    arg_parser = argparse.ArgumentParser(description='');
    #arg_parser.add_argument("coords_input", help="Directory to show-coords tab-delimited input file.");
    arg_parser.add_argument("-g", "--genbank", help="genbank")
    arg_parser.add_argument("-d", "--database", help="database")
    arg_parser.add_argument("-o", "--out_name", help="output name", default="circos_out")
    arg_parser.add_argument("-p", "--phast_database", action="store_true",help="PHAST database")
    arg_parser.add_argument("-e", "--cutoff_evalue", help="e-value cutoff", default=0.00005, type=float)
    arg_parser.add_argument("-c", "--cutoff_coverage", help="query coverage cutoff", default=0.6, type=float)
    arg_parser.add_argument("-i", "--cutoff_identity", help="Percentage identity cutoff", default=50, type=float)
    arg_parser.add_argument("-b", "--blast", help="execute blast", action="store_false")
    arg_parser.add_argument("-s", "--samtools_depth", help="add samtools depth", default=False)
    arg_parser.add_argument("-bt", "--blast_tab", help="uses input blast tab file")
    args = arg_parser.parse_args()

    ##Run main
    if args.phast_database:
        p = PHAST(args.database)
        accession2classif = p.accession2classification
    else:
        accession2classif = False


    print 'execute_blast:', args.blast
    circosf = Genbank2circos(args.genbank,
                             args.database,
                             accession2classification=accession2classif,
                             identity_cutoff=args.cutoff_identity,
                             query_coverage_cutoff=args.cutoff_coverage,
                             evalue_cutoff=args.cutoff_evalue,
                             execute_blast=args.blast,
                             samtools_depth=args.samtools_depth,
                             blast_tab_file=args.blast_tab)

    circosf.write_circos_files(circosf.config, circosf.brewer_conf)
    circosf.run_circos(out_prefix=args.out_name)



