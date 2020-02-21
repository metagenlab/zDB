#! /usr/bin/python

# produce one ptt / record present in the genbank file

from Bio import SeqIO
from optparse import OptionParser
import re
from Bio import SeqUtils

def circos_fasta_draft(fasta_file_name):
    from Bio import SeqIO
    contigs = []
    start = 1
    handle = open(fasta_file_name)
    for seq_record in SeqIO.parse(handle, "fasta"):
        stop = start + len(seq_record.seq)
        contigs.append([seq_record.name, start, stop])
        start = stop+1
    return contigs


def circos_fasta_draft_misc_features(record):
    gap_locations = []
    for feature in record.features:
        if feature.type == "assembly_gap":
            gap_locations.append(feature.location)
    if len(gap_locations) == 0:
        return None
    contigs = []
    for i in range(0, len(gap_locations)):
        if i == 0:
            contigs.append([record.name + "_%s" % 1, 0, int(gap_locations[i].start)])
        else:
            contigs.append([record.name + "_%s" % (i+1), int(gap_locations[i-1].end), int(gap_locations[i].start)])
    contigs.append([record.name + "_%s" % (i+2), int(gap_locations[-1].end), int(len(record.seq))])
    return contigs


def find_index(pattern, seq):
  """Return first item in sequence where f(item) == True."""
  for item in seq:
    if re.match(pattern,item):
      return seq.index(item)

#locus_superantigens = ["saureusC1_00457","saureusC1_00483","saureusC1_00484","saureusC1_00485","saureusC1_00486","saureusC1_00487","saureusC1_00488","saureusC1_00489","saureusC1_00490","saureusC1_00491","saureusC1_00492","saureusC1_00495","saureusC1_01141","saureusC1_01142","saureusC1_01143","saureusC1_01644","saureusC1_02008","saureusC1_02688","saureusC1_02690","saureusC1_02691"]


locus_superantigens = ["CHUV02684",
"CHUV02004",
"CHUV02685",
"CHUV02682",
"CHUV00457",
"CHUV00489",
"CHUV00490",
"CHUV00491",
"CHUV01642",
"CHUV00494",
"CHUV00488",
"CHUV00487",
"CHUV00492",
"CHUV00483",
"CHUV00484",
"CHUV00486",
"CHUV00485",
"CHUV01141",
"CHUV01140",
"CHUV01139"]




"""
class Feature:
    def __init__(self):

        self.contig = "-"
        self.feature_type = "-"
        self.start = "-"
        self.stop = "-"
        self.length = "-"
        self.GC = "-"
        self.sequence = "-"
        self.strand = "-"
        self.gene = "-"
        self.function = "-"
        self.gi = "-"
        self.locus = "-"
        self.protein_id = "-"
        self.inference = "-"
        self.translation = "-"
        self.seq = "-"
        self.id = "-"
        self.product = "-"
        self.orthogroup = "-"



class Record:
    def __init__(self, record):
        self.seq = record.seq
        print record
        #print "length", len(self.seq)
        self.contig =  record.id
        #print "name", self.contig
        self.features = self.get_one_record_features(record)

    def get_one_record_features(self, one_record):

        feature_list = [None] * len(one_record.features)
        for i in range(0,len(one_record.features)):
            #print one_record.features[i]
            if one_record.features[i].type == "misc_feature":
                feature_list[i] = Feature()
                continue
            feature_list[i] = Feature()
            #print  one_record.features[i]

            feature_list[i].feature_type = one_record.features[i].type
            feature_list[i].contig = one_record.id
            feature_list[i].start = one_record.features[i].location.start
            feature_list[i].stop = one_record.features[i].location.end
            feature_list[i].length = len(one_record.features[i].location)
            feature_list[i].strand = one_record.features[i].strand
            try:
                feature_list[i].gene = one_record.features[i].qualifiers['gene'][0]
            except:
                pass
            try:
                gi_position=find_index("GO*",one_record.features[i].qualifiers['db_xref'])
                feature_list[i].gi = one_record.features[i].qualifiers['db_xref'][gi_position][3:]
            except:
                pass

            #geneID= one_record.features[i].qualifiers['db_xref'][1][7:]
            try:
                feature_list[i].locus = one_record.features[i].qualifiers['locus_tag'][0]
            except:
                pass

            try:

                feature_list[i].orthogroup = one_record.features[i].qualifiers['orthogroup'][0]
            except:
                pass
            try:

                feature_list[i].protein_id = one_record.features[i].qualifiers['protein_id'][0]
            except:
                pass
            try:
                feature_list[i].product = one_record.features[i].qualifiers['product'][0]
            except:
                pass
            try:

              feature_list[i].inference = one_record.features[i].qualifiers['inference']
            except:
              pass
            try:
                feature_list[i].translation = one_record.features[i].qualifiers['translation'][0]
            except:
                pass
            feature_list[i].seq = one_record.features[i].extract(self.seq)

            feature_list[i].GC = SeqUtils.GC(feature_list[i].seq)


        return feature_list
"""

def print_circos_record_file(record_list, out_name = "circos_contigs.txt", draft_data = False, draft_coordinates=False):
    i = 1
    f = open(out_name, "w")
    x = 0

    for record in record_list:

      try:
          for contig in draft_data[x]:
              if contig[2] < contig[1]:
                  #print 'impossible contig limits!', contig
                  continue
              #print contig
              if contig[1] == contig[2]:
                  continue
              if draft_coordinates:
                  if i%2 == 0:
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], 0, contig[2]-contig[1], 3))
                  else:
                      ##print "chr -", record.contig, record.contig, "0",len(record.seq), "spectral-5-div-%s" % (4)
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], 0, contig[2]-contig[1], 4))
                  i+=1
              else:
                  if i%2 == 0:
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], contig[1], contig[2], 3))
                  else:
                      ##print "chr -", record.contig, record.contig, "0",len(record.seq), "spectral-5-div-%s" % (4)
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], contig[1], contig[2], 4))
                  i+=1

          x+=1
      except:
            if i%2 == 0:
                f.write("chr - %s %s 0 %s spectral-5-div-%s\n" % (record.id, record.id, len(record.seq), 3))
            else:
                ##print "chr -", record.contig, record.contig, "0",len(record.seq), "spectral-5-div-%s" % (4)
                f.write("chr - %s %s 0 %s spectral-5-div-%s\n" % (record.id, record.id, len(record.seq), 4))
            i += 1
            x += 1
    f.close()




def print_circos_labels_file(record_list, locus2label, out_name = "circos_labels.txt",
                           draft_data=[None], draft_coordinates=False):

    import numpy

    f = open(out_name, "w")

    for y, record in enumerate(record_list):
        for feature in record.features:
            if feature.type == "CDS":

                    try:
                        for i in draft_data[y]:
                            # determine to which contig the feature belong

                            if feature.location.start >= i[1] and feature.location.end <= i[2]:
                                if draft_coordinates:
                                    contig = i[0]
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    contig = i[0]
                                    start = feature.location.start
                                    end = feature.location.end
                    except IndexError:
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end

                    # absurd features
                    if numpy.abs(feature.location.start-feature.location.end) > 50000:
                        continue

                    if 'pseudo' in feature.qualifiers:
                        continue

                    if feature.qualifiers['orthogroup'][0] in locus2label or feature.qualifiers['locus_tag'][0] in locus2label:

                            f.write('%s %s %s %s color=spectral-5-div-2,z=1\n' % (contig,
                                                                               start,
                                                                               end,
                                                                               locus2label[feature.qualifiers['locus_tag'][0]]))

    f.close()



def print_circos_gene_file(record_list, feature_type="CDS", strand ="1",
                           out_name = "circos_genes_plus.txt",
                           locus_highlight=[],
                           draft_data=[None],
                           group_id2orthologs_presence=False,
                           query_taxon_id=False,
                           color_missing=True,
                           draft_coordinates=False,
                           locus_highlight2=[]):

    #print "highlight:", locus_highlight

    import numpy

    if strand == "1":
        f = open(out_name, "w")
    if strand == "-1":
        f = open(out_name, "w")

    #print 'n records:', len(record_list)

    for y, record in enumerate(record_list):
        for feature in record.features:
            if feature.type == feature_type:
                if str(feature.strand) == strand:
                    try:
                        for i in draft_data[y]:
                            # determine to which contig the feature belong

                            if feature.location.start >= i[1] and feature.location.end <= i[2]:
                                if draft_coordinates:
                                    contig = i[0]
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    contig = i[0]
                                    start = feature.location.start
                                    end = feature.location.end
                    except IndexError:
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end

                    '''
                    # in case the second record is not fragmented in several contigs (no draft data)
                    except IndexError:
                        #print "no draft for", record.id
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        #print "no draft for", record.id
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    '''
                    #print "longueur", len(feature.location), type(len(feature.location)), feature.location.start, feature.location.end


                    if numpy.abs(feature.location.start-feature.location.end) > 50000:
                        continue

                    if 'pseudo' in feature.qualifiers:

                        continue
                    try:
                        a = feature.qualifiers['orthogroup']
                        orthology_tag = True
                    except:
                        orthology_tag = False

                    if orthology_tag:
                        if feature.qualifiers['orthogroup'][0] in locus_highlight or feature.qualifiers['locus_tag'][0] in locus_highlight:

                            #print "COLORS!!!!!!!!!!!!!!!!!!!"
                            try:
                                '''
                                f.write('%s %s %s fill_color=spectral-5-div-4,id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                                         start,
                                                                                                         end,
                                                                                                         feature.qualifiers['locus_tag'][0],
                                                                                                         feature.qualifiers['gene'][0],
                                                                                                         re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                                         feature.qualifiers['product'][0])))
                                '''
                                f.write('%s %s %s fill_color=piyg-5-div-1,id=%s,z=1\n' % (contig,
                                                                        start,
                                                                        end, feature.qualifiers['locus_tag'][0]))


                            except:
                                '''
                                f.write('%s %s %s fill_color=spectral-5-div-4,id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                                         start,
                                                                                                         end,
                                                                                                         feature.qualifiers['locus_tag'][0],
                                                                                                         "-",
                                                                                                         re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                                         feature.qualifiers['product'][0])))
                                '''
                                f.write('%s %s %s fill_color=spectral-5-div-4,id=%s,z=1\n' % (contig,
                                                                        start,
                                                                        end,
                                                                        feature.qualifiers['locus_tag'][0]))
                        elif feature.qualifiers['orthogroup'][0] in locus_highlight2 or feature.qualifiers['locus_tag'][0] in locus_highlight2:
                            try:
                                f.write('%s %s %s fill_color=paired-11-qual-10,id=%s,z=1\n' % (contig,
                                                                        start,
                                                                        end, feature.qualifiers['locus_tag'][0]))


                            except:
                                f.write('%s %s %s fill_color=spectral-5-div-4,id=%s,z=1\n' % (contig,
                                                                        start,
                                                                        end,
                                                                        feature.qualifiers['locus_tag'][0]))


                        else:

                            f.write('%s %s %s id=%s\n' % (contig,
                                                  start,
                                                  end, feature.qualifiers['locus_tag'][0]))

                    elif feature.qualifiers['locus_tag'][0] in locus_highlight:
                        try:

                            f.write('%s %s %s fill_color=piyg-5-div-1,id=%s,z=1\n' % (contig,
                                                                    start,
                                                                    end, feature.qualifiers['locus_tag'][0]))


                        except:
                            f.write('%s %s %s fill_color=spectral-5-div-4,id=%s,z=1\n' % (contig,
                                                                    start,
                                                                    end,
                                                                    feature.qualifiers['locus_tag'][0]))

                        else:

                            f.write('%s %s %s id=%s\n' % (contig,
                                                  start,
                                                  end, feature.qualifiers['locus_tag'][0]))
                    # not in locus_highlight
                    else:
                        if group_id2orthologs_presence and query_taxon_id and color_missing:
                            try:
                                pseudo = feature.qualifiers['pseudogene']
                                #print 'pseudogene, continue'
                                continue
                            except:
                                pass
                            if group_id2orthologs_presence[feature.qualifiers['orthogroup'][0]][int(query_taxon_id)] == 0:
                                f.write('%s %s %s fill_color=orange\n' % (contig,
                                                                        start,
                                                                        end))
                            else:
                                f.write('%s %s %s id=%s\n' % (contig,
                                                      start,
                                                      end, feature.qualifiers['locus_tag'][0]))
                        else:

                            f.write('%s %s %s id=%s\n' % (contig,
                                                  start,
                                                  end, feature.qualifiers['locus_tag'][0]))

                        '''
                            try:
                                f.write('%s %s %s fill_color=orange, id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                                             start,
                                                                                                             end,
                                                                                                             feature.qualifiers['locus_tag'][0],
                                                                                                             feature.qualifiers['gene'][0],
                                                                                                             re.sub("[ |\-|(|)|\]|\[|\.|,]+", "_",
                                                                                                                    feature.qualifiers['product'][0])))
                            except:
                                f.write('%s %s %s fill_color=orange, id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                                             start,
                                                                                                             end,
                                                                                                             feature.qualifiers['locus_tag'][0],
                                                                                                             "-",
                                                                                                             re.sub("[ |\-|(|)|\]|\[|\.|,]+", "_",
                                                                                                                    feature.qualifiers['product'][0])))

                        else:

                            try:
                                f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                     start,
                                                                                     end,
                                                                                     feature.qualifiers['locus_tag'][0],
                                                                                     feature.qualifiers['gene'][0],
                                                                                     re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                            feature.qualifiers['product'][0])))
                            except:
                                f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                     start,
                                                                                     end,
                                                                                     feature.qualifiers['locus_tag'][0],
                                                                                     '-',
                                                                                     re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                            feature.qualifiers['product'][0])))
                    else:

                        try:
                            f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                 start,
                                                                                 end,
                                                                                 feature.qualifiers['locus_tag'][0],
                                                                                 feature.qualifiers['gene'][0],
                                                                                 re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                        feature.qualifiers['product'][0])))
                        except:
                            try:
                                f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                 start,
                                                                                 end,
                                                                                 feature.qualifiers['locus_tag'][0],
                                                                                 '-',
                                                                                 re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                        feature.qualifiers['product'][0])))
                            except:
                                f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                 start,
                                                                                 end,
                                                                                 feature.qualifiers['locus_tag'][0],
                                                                                 '-',
                                                                                 re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                        "-")))
                        '''

            if feature.type == 'rRNA':
                #print 'rrna--------------------'
                if str(feature.strand) == strand:
                    try:
                        for i in draft_data[y]:
                            # determine to which contig the feature belong

                            if feature.location.start >= i[1] and feature.location.end <= i[2]:
                                if draft_coordinates:
                                    contig = i[0]
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    contig = i[0]
                                    start = feature.location.start
                                    end = feature.location.end
                    except IndexError:
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end

                    #print 'rrna position:', contig, start, end, feature.qualifiers
                    f.write('%s %s %s fill_color=pblue\n ' % (contig,
                                                            start,
                                                            end))


            if feature.type == 'tRNA':
                if str(feature.strand) == strand:
                    try:
                        for i in draft_data[y]:
                            # determine to which contig the feature belong

                            if feature.location.start >= i[1] and feature.location.end <= i[2]:
                                if draft_coordinates:
                                    contig = i[0]
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    contig = i[0]
                                    start = feature.location.start
                                    end = feature.location.end
                    except IndexError:
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end

                    f.write('%s %s %s fill_color=pred\n ' % (contig,
                                                            start,
                                                            end))


    f.close()
"""
                else :
                    f.write('%s %s %s id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                          feature.location.start,
                                                                          feature.location.end,
                                                                          feature.qualifiers['locus_tag'][0],
                                                                          feature.qualifiers['gene'][0],
                                                                          re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                 feature.qualifiers['product'][0])))

"""




def print_circos_GC_file(record_list, feature_type="CDS", out_directory=""):


    """
    print "draft_data", draft_data
    f = open(out_name, "w")
    for record in record_list:
    for feature in record.features:
      if feature.feature_type == feature_type:
        if draft_data:
            for i in draft_data:
              if feature.start >= i[1] and feature.stop <= i[2]:
                contig = i[0]
        else:
            contig = feature.contig
        f.write("%s %s %s %s\n" % (contig, feature.start, feature.stop, feature.GC))
    """
    from chlamdb.plots import GC
    import os

    out_var_file = os.path.join(out_directory, 'circos_GC_var_%s.txt' % record_list[0].name)
    out_skew_file = os.path.join(out_directory, 'circos_GC_skew_%s.txt' % record_list[0].name)
    f = open(out_var_file, 'w')
    g = open(out_skew_file, 'w')

    out_var = ''
    out_skew = ''
    for record in record_list:
        # this function handle scaffolds (split sequence when encountering NNNNN regions)
        out_var += GC.circos_gc_var(record, windows=1000)
        out_skew += GC.circos_gc_skew(record, windows=1000)
    f.write(out_var)
    g.write(out_skew)



    return (out_var_file, out_skew_file)


def print_blasnr_circos_files(record_list, db_name, out_directory, draft_coordinates=False, exclude_family=False, taxon_list =False):

    '''
    :param record:
    :return: circos string with difference as compared to the average GC
    ex: average = 32
        GC(seq[3000:4000]) = 44
        diff = 44 - 32 = 12%
    if taxon_list: restrict counts of homologs in the set if specified taxons

    '''
    import numpy
    import pandas
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    #print '##########  %s  ###############' % exclude_family

    accessions = []
    for record in record_list:
        accessions.append(record.id.split('.')[0])

    from chlamdb.biosqldb import biosql_own_sql_tables as bsql

    #print '-------------- taxon list ------------'
    #print taxon_list

    if not taxon_list:
        locus_tag2n_genomes_dico = bsql.locus_tag2presence_in_n_genomes(db_name)
    else:
        taxon_list = [str(i) for i in taxon_list]
        filter = '`' + '`,`'.join(taxon_list) + '`'
        sql = 'select orthogroup,%s from comparative_tables_orthology' % (filter)

        data = numpy.array([list(i) for i in server.adaptor.execute_and_fetchall(sql,)])
        count_df = pandas.DataFrame(data, columns=['orthogroup'] + taxon_list)
        count_df = count_df.set_index(['orthogroup'])
        count_df = count_df.apply(pandas.to_numeric, args=('coerce',))

        group2n_genomes = (count_df > 0).sum(axis=1)
        #print 'group2n_genomes--------------'
        #print group2n_genomes

        accession_filter = '"'+'","'.join(accessions) + '"'
        sql = 'select locus_tag, orthogroup from orthology_detail where accession in (%s)' % (accession_filter)
        locus2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        locus_tag2n_genomes_dico = {}
        for locus in locus2orthogroup:
            locus_tag2n_genomes_dico[locus] = group2n_genomes[locus2orthogroup[locus]]

        #locus_tag2n_genomes = manipulate_biosqldb.get_orthology_table_subset()

    locus_tag2n_blastnr_dico = {}
    locus_tag2n_blast_bacteria = {}
    locus_tag2n_blast_eukariota = {}
    locus_tag2n_archae = {}
    locus_tag2n_chlamydiae = {}
    locus_tag2n_non_chlamydiae = {}
    locus_tag2n_paralogs = {}

    #print 'getting dictionnaries'
    for accession in accessions:

        #print "locus_tag2n_nr_hits"
        locus_tag2n_blastnr_dico.update(bsql.locus_tag2n_nr_hits(db_name, accession, exclude_family=exclude_family))
        #print "locus_tag2n_blast_superkingdom bacteria"
        locus_tag2n_blast_bacteria.update(bsql.locus_tag2n_blast_superkingdom(db_name, accession, superkingdom="Bacteria", exclude_family=exclude_family))
        #print "locus_tag2n_blast_superkingdom euk"
        locus_tag2n_blast_eukariota.update(bsql.locus_tag2n_blast_superkingdom(db_name, accession, superkingdom="Eukaryota"))
        #print "locus_tag2n_blast_superkingdom archae"
        locus_tag2n_archae.update(bsql.locus_tag2n_blast_superkingdom(db_name, accession, superkingdom="Archaea"))
        #print "locus_tag2n_blast_bacterial_phylum chlamydiae reverse = F"
        locus_tag2n_chlamydiae.update(bsql.locus_tag2n_blast_bacterial_phylum(db_name, accession, phylum="Chlamydiae", reverse=False, exclude_family=exclude_family))
        #print "locus_tag2n_blast_bacterial_phylum chlamydiae reverse = T"
        locus_tag2n_non_chlamydiae.update(bsql.locus_tag2n_blast_bacterial_phylum(db_name, accession, phylum="Chlamydiae", reverse=True, exclude_family=exclude_family))
        #print "locus_tag2n_paralogs"
        locus_tag2n_paralogs.update(bsql.locus_tag2n_paralogs(db_name, accession))


    circos_string_n_genome_presence = ''
    circos_string_n_blastnr = ''
    circos_string_n_blast_bacteria = ''
    circos_string_n_blast_eukariota = ''
    circos_string_n_archae = ''
    circos_string_n_chlamydiae = ''
    circos_string_n_non_chlamydiae = ''
    circos_string_n_paralogs = ''
    circos_string_stacked_chlamydiales = ''
    circos_string_stacked_chlamydiales_and_non_bacteria = ''

    draft_data = []
    for biorecord in record_list:
        temp = circos_fasta_draft_misc_features(biorecord)
        draft_data.append(temp)



    y = 0
    for record in record_list:
        for feature in record.features:
            if feature.type == 'CDS':
                if numpy.abs(feature.location.start-feature.location.end) > 50000:
                    continue
                if not 'pseudo' in feature.qualifiers:
                    try:
                        for i in draft_data[y]:
                            # find in which contig the feature is located
                            if int(feature.location.start) >= i[1] and int(feature.location.end) <= i[2]:
                                contig = i[0]
                                # keep contig coordinates or not
                                if draft_coordinates:
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    start = feature.location.start
                                    end = feature.location.end
                    # in case we have draft data for one record and not the other (i.e not for a plasmid)

                    except:
                        contig = record.id
                        start = feature.location.start
                        end = feature.location.end
                    #print "contig", contig

                    try:
                        n_genomes = locus_tag2n_genomes_dico[feature.qualifiers['locus_tag'][0]]
                    except KeyError:
                        import sys
                        #print 'problem with n homologs in', feature
                        sys.exit()

                    try:
                        n_blastnr = locus_tag2n_blastnr_dico[feature.qualifiers['locus_tag'][0]]
                        if n_blastnr == None:
                            n_blastnr = 0
                    except:
                        n_blastnr = 0
                    try:
                        n_blast_bacteria = locus_tag2n_blast_bacteria[feature.qualifiers['locus_tag'][0]]
                    except KeyError:
                        n_blast_bacteria = 0

                    try:
                        n_blast_eukariota = locus_tag2n_blast_eukariota[feature.qualifiers['locus_tag'][0]]
                    except KeyError:
                        n_blast_eukariota = 0
                    try:
                        n_archae = locus_tag2n_archae[feature.qualifiers['locus_tag'][0]]
                    except KeyError:
                        n_archae = 0

                    try:
                        n_chlamydiae = locus_tag2n_chlamydiae[feature.qualifiers['locus_tag'][0]]
                    except KeyError:
                        n_chlamydiae = 0

                    try:
                        n_non_chlamydiae = locus_tag2n_non_chlamydiae[feature.qualifiers['locus_tag'][0]]
                    except KeyError:
                        n_non_chlamydiae = 0

                    try:
                        n_paralogs = locus_tag2n_paralogs[feature.qualifiers['locus_tag'][0]]
                    except KeyError:
                        n_paralogs = 0

                    circos_string_n_genome_presence += "%s %s %s %s id=%s\n" % (contig, start, end, n_genomes,feature.qualifiers['locus_tag'][0])
                    circos_string_n_blastnr += "%s %s %s %s id=%s\n" % (contig, start, end, n_blastnr,feature.qualifiers['locus_tag'][0])
                    circos_string_n_blast_bacteria += "%s %s %s %s id=%s\n" % (contig, start, end, n_blast_bacteria,feature.qualifiers['locus_tag'][0])
                    circos_string_n_blast_eukariota += "%s %s %s %s id=%s\n" % (contig, start, end, n_blast_eukariota,feature.qualifiers['locus_tag'][0])
                    circos_string_n_archae += "%s %s %s %s id=%s\n" % (contig, start, end, n_archae,feature.qualifiers['locus_tag'][0])
                    circos_string_n_chlamydiae += "%s %s %s %s id=%s\n" % (contig, start, end, n_chlamydiae,feature.qualifiers['locus_tag'][0])
                    circos_string_n_non_chlamydiae += "%s %s %s %s id=%s\n" % (contig, start, end, n_non_chlamydiae,feature.qualifiers['locus_tag'][0])
                    circos_string_n_paralogs += "%s %s %s %s id=%s\n" % (contig, start, end, n_paralogs,feature.qualifiers['locus_tag'][0])
                    circos_string_stacked_chlamydiales += "%s %s %s %s,%s id=%s\n" % (contig, start, end, n_chlamydiae, n_non_chlamydiae,feature.qualifiers['locus_tag'][0])
                    circos_string_stacked_chlamydiales_and_non_bacteria += "%s %s %s %s,%s,%s id=%s\n" % (contig, start, end, n_chlamydiae, n_non_chlamydiae,str(int(n_blast_eukariota)+int(n_archae)),feature.qualifiers['locus_tag'][0])
        y += 1
    import os
    all_file_names = {}
    all_file_names['file_n_genomes'] = os.path.join(out_directory,"circos_n_genome_presence.txt")
    all_file_names['file_n_blastnr'] = os.path.join(out_directory,"circos_n_blastnr.txt")
    all_file_names['file_n_blast_bacteria'] = os.path.join(out_directory,"circos_n_blast_bactera.txt")
    all_file_names['file_n_blast_eukaryote'] = os.path.join(out_directory,"circos_n_blast_eukaryota.txt")
    all_file_names['file_n_blast_archae'] = os.path.join(out_directory,"circos_n_blast_archae.txt")
    all_file_names['file_n_blast_chlamydiae'] = os.path.join(out_directory,"circos_n_blast_chlamydiae.txt")
    all_file_names['file_n_blast_non_chlamydiae'] = os.path.join(out_directory,"circos_n_blast_non_chlamydiae.txt")
    all_file_names['file_n_paralogs'] = os.path.join(out_directory,"circos_n_paralogs.txt")
    all_file_names['file_stacked_chlamydiales'] = os.path.join(out_directory,"circos_stacked_chlamydiales.txt")
    all_file_names['file_stacked_chlamydiales_non_prokaryotes'] = os.path.join(out_directory,"circos_stacked_chlamydiales_non_euk.txt")

    all_file_names['gc_var_file'], all_file_names['gc_skew_file'] = print_circos_GC_file(record_list, feature_type="CDS", out_directory=out_directory)

    with open(all_file_names['file_n_genomes'], "w") as f:
        f.write(circos_string_n_genome_presence)

    with open(all_file_names['file_n_blastnr'], "w") as f:
        f.write(circos_string_n_blastnr)

    with open(all_file_names['file_n_blast_bacteria'], "w") as f:
        f.write(circos_string_n_blast_bacteria)

    with open(all_file_names['file_n_blast_eukaryote'], "w") as f:
        f.write(circos_string_n_blast_eukariota)

    with open(all_file_names['file_n_blast_archae'], "w") as f:
        f.write(circos_string_n_archae)

    with open(all_file_names['file_n_blast_chlamydiae'], "w") as f:
        f.write(circos_string_n_chlamydiae)

    with open(all_file_names['file_n_blast_non_chlamydiae'], "w") as f:
        f.write(circos_string_n_non_chlamydiae)

    with open(all_file_names['file_n_paralogs'], "w") as f:
        f.write(circos_string_n_paralogs)

    with open(all_file_names['file_stacked_chlamydiales'], "w") as f:
        f.write(circos_string_stacked_chlamydiales)
    with open(all_file_names['file_stacked_chlamydiales_non_prokaryotes'], "w") as f:
        f.write(circos_string_stacked_chlamydiales_and_non_bacteria)
    #all_file_names.update(gc_files)

    # return blastnr + GC files names
    return all_file_names






def orthology_circos_files(server,
                           record_list,
                           reference_taxon_id,
                           biodatabase_name, out_dir,
                           locus_highlight=[],
                           feature_type="CDS",
                           taxon_list=False,
                           draft_data=[None],
                           query_taxon_id=False,
                           draft_coordinates=False,
                           color_missing=True,
                           locus2label=False,
                           show_homologs=True,
                           get_orthogroup_counts=False,
                           locus_highlight2=[]):

    from chlamdb.biosqldb import biosql_own_sql_tables
    import os
    from chlamdb.biosqldb import orthogroup_identity_db
    from chlamdb.biosqldb import manipulate_biosqldb

    #print "orthology_circos_files"
    #print "draft_data", draft_data
    #print "locus highlight", locus_highlight
    all_file_names = {}
    #ortho_size = mysqldb_load_mcl_output.get_all_orthogroup_size(server, biodatabase_name)

    sql = 'select orthogroup, count(*) from orthology_detail group by orthogroup'
    ortho_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    #print 'getting each orthogroup size'
    #try:
    #    ortho_identity = orthogroup_identity_db.orthogroup2average_identity(biodatabase_name)
    #except:
    #    pass
    #ortho_family_size = mysqldb_load_mcl_output.get_family_size(server, biodatabase_name)

    ortho_table = "orthology_%s" % biodatabase_name
    if get_orthogroup_counts:
        group_id2orthologs_presence = biosql_own_sql_tables.get_orthology_matrix_merging_plasmids(server, biodatabase_name, taxon_list)
    else:
        group_id2orthologs_presence = False
    # get taxons ids used as column names in orthology table
    if taxon_list is False:
        taxon_list = manipulate_biosqldb.get_column_names(server, ortho_table)[1:]

    reference_taxon_id = str(reference_taxon_id)

    all_file_names["contigs"] = os.path.join(out_dir,
                                           "circos_contigs_%s.txt" % reference_taxon_id)
    all_file_names["minus"] = os.path.join(out_dir,
                                         "circos_genes_minus_%s.txt" % reference_taxon_id)
    all_file_names["plus"] = os.path.join(out_dir,
                                        "circos_genes_plus_%s.txt" % reference_taxon_id)
    all_file_names["GC"] = os.path.join(out_dir,
                                      "circos_GC_%s.txt" % reference_taxon_id)
    all_file_names["orthogroups"] = os.path.join(out_dir,
                                               "circos_orthogroup_size_%s.txt"  % reference_taxon_id)
    all_file_names["orthogroups_identity"] = os.path.join(out_dir,
                                               "circos_orthogroup_identity_%s.txt"  % reference_taxon_id)
    all_file_names["orthogroups_family"] = os.path.join(out_dir,
                                               "circos_orthogroup_family_size_%s.txt"  % reference_taxon_id)

    if locus2label:
        all_file_names["labels"] = os.path.join(out_dir,
                                                   "circos_labels_%s.txt"  % reference_taxon_id)

        print_circos_labels_file(record_list,
                                 locus2label,
                                 all_file_names["labels"],
                                 draft_data=draft_data,
                                 draft_coordinates=draft_coordinates)

    #print "writing record file"
    print_circos_record_file(record_list, out_name = all_file_names["contigs"], draft_data=draft_data)
    #print "writing minus strand file"
    print_circos_gene_file(record_list, strand="-1",
                           out_name = all_file_names["minus"],
                           draft_data=draft_data,
                           locus_highlight=locus_highlight,
                           group_id2orthologs_presence=group_id2orthologs_presence,
                           query_taxon_id=query_taxon_id,
                           draft_coordinates=draft_coordinates,
                           locus_highlight2=locus_highlight2)

    #print "writing plus strand file"
    print_circos_gene_file(record_list,
                           strand="1",
                           out_name = all_file_names["plus"],
                           draft_data=draft_data,
                         locus_highlight=locus_highlight,
                           group_id2orthologs_presence=group_id2orthologs_presence,
                         query_taxon_id=query_taxon_id, draft_coordinates=draft_coordinates,
                           locus_highlight2=locus_highlight)
    #print "writing GC file"



    all_file_names["GC_var"], all_file_names["GC_skew"] = print_circos_GC_file(record_list, out_directory = out_dir)

    locus_highlight = []

    f = open(all_file_names["orthogroups"], "w")
    #g = open(all_file_names["orthogroups_identity"], "w")
    #h = open(all_file_names["orthogroups_family"], "w")

    taxon_files = []
    all_file_names["genomes"] = []

    temp_taxon_list = [int(i) for i in taxon_list]
    for i in range(0, len(taxon_list)):
        #print i
        if int(taxon_list[i]) == int(reference_taxon_id):
            temp_taxon_list.pop(temp_taxon_list.index(int(reference_taxon_id)))
        else:
            file_name = os.path.join(out_dir, "circos_taxon_%s_vs_%s.txt" % (reference_taxon_id, taxon_list[i]))
            all_file_names["genomes"].append(file_name)
            taxon_files.append(open(file_name, "w"))
    #print 'get closest homolog identity dictionnary'
    if show_homologs:
        locus2locus_identity = biosql_own_sql_tables.circos_locus2taxon_highest_identity(biodatabase_name,
                                                                                         reference_taxon_id,
                                                                                         use_identity_closest_homolog2_table=True)

    #print "get plasmid locus"
    plasmid_locus = biosql_own_sql_tables.get_locus2plasmid_or_not(biodatabase_name)

    taxon_list = temp_taxon_list
    y = 0
    for record in record_list:
        for feature in record.features:
            if feature.type == feature_type:
                if not 'pseudo' in feature.qualifiers:
                    try:
                        for i in draft_data[y]:
                            # cas des locations splittes en 2
                            if feature.location.__class__.__name__ != "CompoundLocation":
                                if int(feature.location.start) >= i[1] and int(feature.location.end) <= i[2]:
                                    #print "OKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK"
                                    contig = i[0]
                                    if draft_coordinates:
                                        start = feature.location.start - i[1]
                                        end = feature.location.end - i[1]
                                    else:
                                        start = feature.location.start
                                        end = feature.location.end
                            else:
                                if draft_coordinates:
                                    start = feature.location.start - i[1]
                                    end = feature.location.parts[0].end - i[1]
                                else:
                                    start = feature.location.start
                                    end = feature.location.parts[0].end

                    except:
                        # cas des locations splittes en 2
                        if feature.location.__class__.__name__ != "CompoundLocation":
                            contig = record.id
                            start = feature.location.start
                            end = feature.location.end
                        else:
                            contig = record.id
                            start = feature.location.start
                            end = feature.location.parts[0].end
                    #print "contig", contig

                    try:
                        line = "%s %s %s %s\n" % (contig, start, end, ortho_size[feature.qualifiers['orthogroup'][0]])
                    except:
                        #print "problem ith", feature
                        #print contig, start, end
                        #print feature
                        continue
                    #try:
                    #    line2 = "%s %s %s %s\n" % (contig, start, end, ortho_identity[feature.qualifiers['orthogroup'][0]][0])
                    #except:
                    #    pass
                    #line3 = "%s %s %s %s\n" % (contig, start, end, ortho_family_size[feature.qualifiers['orthogroup'][0]])
                    f.write(line)
                    #try:
                    #    g.write(line2)
                    #except UnboundLocalError:
                    #    pass
                    #h.write(line3)
                    #print "################## highlight ###################"
                    #print locus_highlight


                    for taxon_id, taxon_file in zip(taxon_list, taxon_files):
                        #print group_id2orthologs_presence[feature.qualifiers['orthogroup'][0]]
                        if group_id2orthologs_presence[feature.qualifiers['orthogroup'][0]][int(taxon_id)] > 0:
                            #if locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][1] == 'PMK1_b00074':
                            #    print plasmid_locus[locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][1]]
                            #    print plasmid_locus[locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][1]] == 'True'
                            if plasmid_locus[locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][1]] != 'True':
                                if feature.qualifiers['orthogroup'][0] not in locus_highlight and feature.qualifiers['locus_tag'][0] not in locus_highlight:

                                    taxon_file.write("%s %s %s %s id=%s\n" % (contig,
                                                                              start,
                                                                              end,
                                                                              locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][0],
                                                                              locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][1]))

                                else:
                                    #print "COLORS!!!!!!!!!!!!!!!!!!! ----- taxon" # pink piyg-5-div-1 spectral-5-div-4 green
                                    taxon_file.write("%s %s %s %s color=piyg-5-div-1,id=%s,z=2\n" % (contig,
                                                                                                 start,
                                                                                                 end,
                                                                                                 locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][0],
                                                                                                 locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][1]))
                            else:
                                if feature.qualifiers['orthogroup'][0] not in locus_highlight and feature.qualifiers['locus_tag'][0] not in locus_highlight:

                                    taxon_file.write("%s %s %s %s id=%s,color=accent-4-qual-4\n" % (contig,
                                                                              start,
                                                                              end,
                                                                              locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][0],
                                                                              locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][1]))

                                else:
                                    #print "COLORS!!!!!!!!!!!!!!!!!!! ----- taxon"
                                    taxon_file.write("%s %s %s %s color=piyg-5-div-1,id=%s,color=accent-4-qual-4,z=2\n" % (contig,
                                                                                                 start,
                                                                                                 end,
                                                                                                 locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][0],
                                                                                                 locus2locus_identity[feature.qualifiers['locus_tag'][0]][taxon_id][1]))




                        elif color_missing:
                            #if feature.qualifiers['orthogroup'][0] not in locus_highlight:
                            #    print 'no colors'
                            taxon_file.write("%s %s %s %s\n" % (contig, start, end, 0))
                            #else:
                            #    print "COLORS!!!!!!!!!!!!!!!!!!! ----- taxon"
                            #    taxon_file.write("%s %s %s %s fill_color=spectral-5-div-4,z=2\n" % (contig, start, end, 0))
                        else:
                            pass #taxon_file.write("%s %s %s\n" % (contig, feature.start, feature.stop))
        y += 1
    return all_file_names

  #print manipulate_biosqldb.get_taxon_id_list(server, biodatabase_name)


class Circos_config:
  def __init__(self, caryotype_file,
               chr_spacing_list=[],
               show_ticks="yes",
               show_tick_labels="yes",
               ideogram_spacing=0,
               radius=0.75,
               label_radius=0.05,
               show_ideogram_labels="no",
               color_files=""):
    import re

    self.plots = ""
    self.links = ""
    self.highlights = ""

    self.template_caryotype= "karyotype = %s\n" \
                             " chromosomes_units           = 10000\n" \
                             " chromosomes_display_default = yes\n" % caryotype_file


    self.template_ideograms = "<ideogram>\n" \
                              " <spacing>\n" \
                              " default            = %su\n" \
                              " %s" \
                              " </spacing>\n" \
                              " \n" \
                              " # thickness and color of ideograms\n" \
                              " thickness          = 12p\n" \
                              " stroke_thickness   = 1\n" \
                              " stroke_color       = black\n" \
                              " \n" \
                              " # the default chromosome color is set here and any value\n" \
                              " # defined in the karyotype file overrides it\n" \
                              " fill               = yes\n" \
                              " fill_color         = black\n" \
                              " \n" \
                              " # fractional radius position of chromosome ideogram within image\n" \
                              " radius             = %sr\n" \
                              " show_label         = %s\n" \
                              " label_font         = default\n" \
                              " label_radius       = dims(ideogram,radius) + %sr\n" \
                              " label_size         = 30\n" \
                              " label_parallel     = no\n" \
                              " \n" \
                              " # show_bands determines whether the outline of cytogenetic bands\n" \
                              " # will be seen\n" \
                              " show_bands         = yes\n" \
                              " band_stroke_thickness = 1\n" \
                              " " \
                              " # in order to fill the bands with the color defined in the karyotype\n" \
                              " # file you must set fill_bands\n" \
                              " fill_bands         = yes\n" \
                              " band_transparency  = 1\n" \
                              " \n" \
                              " </ideogram>\n" % (ideogram_spacing,
                                                  self.add_spacing(chr_spacing_list),
                                                  radius,
                                                  show_ideogram_labels,
                                                  label_radius)

    self.end_ticks = "  <tick>\n" \
                      "  multiplier   = 1\n" \
                      "  position = end\n" \
                      " show_ticks         = yes\n" \
                      "  size              = 4p\n" \
                      "  show_label        = no\n" \
                      "  label_size        = 0p\n" \
                      "  format    = %s bp\n" \
                      "  </tick>\n" \

    self.big_ticks = " <tick>\n" \
                          " show_ticks         = yes\n" \
                          " skip_first_label = no\n" \
                          " multiplier   = 10/1u\n" \
                          " spacing           = 10u\n" \
                          " size              = 15p\n" \
                          " show_label        = yes\n" \
                          " label_size        = 25p\n" \
                          " format            = %s kb\n" \
                          " thickness         = 2p\n" \
                          " </tick>\n" \


    self.template_ticks = "show_ticks         = %s\n" \
                          " show_tick_labels   = %s\n" \
                          " \n" \
                          " <ticks>\n" \
                          " tick_label_font    = condensed\n" \
                          " radius             = dims(ideogram,radius_outer)\n" \
                          " label_offset       = 8p\n" \
                          " label_size         = 4p\n" \
                          " color              = black\n" \
                          " thickness          = 2p\n" \
                          " \n" \
                          " \n" \
                          " <tick>\n" \
                          " show_ticks         = yes\n" \
                          " skip_first_label = no\n" \
                          " multiplier   = 10/1u\n" \
                          " spacing           = 1u\n" \
                          " size              = 5p\n" \
                          " show_label        = no\n" \
                          " label_size        = 5p\n" \
                          " format            = %s\n" \
                          " </tick>\n" \
                          " %s\n" \
                          " \n" \
                          " \n" \
                          " </ticks>\n" % (show_ticks, show_tick_labels, "%%.1d", self.big_ticks % "%%.1d")
    print (self.template_ticks)




    self.template_rules = "<rules>\n" \
                          " %s\n" \
                          "</rules>\n"

    self.template_backgrounds = "<backgrounds>\n" \
                          " %s\n" \
                          "</backgrounds>\n"

    # #" <<include colors.rn.conf>>\n" \
    self.settings ="<colors>\n" \
                   " %s\n" \
                   " #<<include brewer.all.conf>>\n" \
                   " </colors>\n" \
                   " <image>\n"\
                   " image_map_use      = yes\n" \
                   " image_map_overlay  = no\n" \
                   " image_map_overlay_stroke_color     = red\n" \
                   " <<include image.conf>>\n" \
                   " </image>\n" \
                   " #includes  etc/colors.conf\n" \
                   " #          etc/fonts.conf\n" \
                   " #          etc/patterns.conf\n" \
                   " <<include colors_fonts_patterns.conf>>\n" \
                   " # system and debug settings\n" \
                   " <<include housekeeping.conf>>\n" \
                   " anti_aliasing*     = no\n"

    self.complete_file = self.template_caryotype + self.template_ideograms + self.template_ticks + "%s %s %s" + self.settings % (color_files)

  def _template_spacing(self, chr1, chr2):
    template = '<pairwise %s %s>\n' \
               ' spacing = 2u\n' \
               '</pairwise>\n' % (chr1, chr2)
    return template

  def _template_plot(self, file, type="line", r0=1,
                     r1=1.05, color=False, fill_color="red", thickness = "0.8p", z = 1, rules ="", backgrounds="", url="", min=False, max=False):

        #print 'template color------------------', color
        template1 = "<plot>\n" \
               "type		    = %s\n" \
               " url                = %s[id]\n" \
               " r0                 = %s\n" \
               " r1                 = %s\n" \
               " fill_color         = %s\n" \
               " thickness          = %s\n" \
               " file               = %s\n" \
               " z                  = %s\n" % (type, url, r0, r1, fill_color, thickness, file, z)
        #print '--------------- COLOR--------------'
        #print color
        if color:
            template1+= " color          = %s\n" % color
        if min:
            template1+= " min          = %s\n" % min
        if max:
            max+= " max          = %s\n" % max
        template_rules = " %s\n" \
               " %s\n" \
               " </plot>\n" % (rules, backgrounds)

        return template1 + template_rules

  def template_rule(self, condition, fill_color):
    template = "<rule>\n" \
               " condition          = %s\n" \
               " fill_color         = %s\n" \
               " </rule>\n" % (condition, fill_color)
    return template

  def template_background(self, color):
    template = "<background>\n" \
               " color         = %s\n" \
               " </background>\n" % (color)
    return template



  def _template_link(self, link_file, color="black_a5", thickness=1):
    template = "<link>\n" \
               "ribbon = yes\n" \
               "file          = %s\n" \
               "color         = %s\n" \
               "radius        = 0.99r\n" \
               "bezier_radius = 0.1r\n" \
               "thickness     = %s\n" \
               " </link>\n" % (link_file, color, thickness)
    return template


  def _template_highlight(self, file, fill_color="grey_a1", r1="1.55r", r0="1.50r", url="/chlamdb/locusx/chlamydia_03_15/"):
    template ="<highlight>\n" \
              " fill_color = %s\n" \
              " file       = %s\n" \
              " r1         = %s\n" \
              " r0         = %s\n" \
              " url = %s[id]\n" \
              " </highlight>\n" % (fill_color, file, r1, r0, url)

    return template


  def add_plot(self, file, type="line", r0=1, r1=1.05,
               fill_color="grey_a1", thickness = "2p", z = 1, rules ="", backgrounds="", url="", min=False, max=False, color=False):
    plot = self._template_plot(file, type, r0, r1, color, fill_color, thickness, z, rules, backgrounds, url, min=min, max=max)
    if len(re.findall("</plots>", self.plots))>0:
      # remove end balise
      self.plots = re.sub("</plots>", "", self.plots)
      # add new plot and end balise
      self.plots = self.plots + plot + "</plots>\n"
    else:
      self.plots = "<plots>\n" + plot + "</plots>\n"



  def add_link(self, link_file, color="black_a5", thickness=1):
    link = self._template_link(link_file, color=color, thickness=thickness)
    if len(re.findall("</links>", self.plots))>0:
      # remove end balise
      self.links = re.sub("</links>", "", self.plots)
      # add new plot and end balise
      self.links = self.links + link + "</links>\n"
    else:
      self.links = "<links>\n" + link + "</links>\n"


  def add_highlight(self, file, fill_color="grey_a1", r1="1.55r", r0="1.50r", href=""):
    highlight = self._template_highlight(file, fill_color, r1, r0, href)
    if len(re.findall("</highlights>", self.highlights))>0:
      # remove end balise
      self.highlights = re.sub("</highlights>", "", self.highlights)
      # add new plot and end balise
      self.highlights = self.highlights + highlight + "</highlights>\n"
    else:
      self.highlights = "<highlights>\n" + highlight + "</highlights>\n"

  def add_spacing(self, chr_pair_list):
      all_spacing = ''
      if len(chr_pair_list) == 0:
          return ''
      else:
          for pair in chr_pair_list:
              all_spacing += self._template_spacing(pair[0], pair[1])
          return all_spacing

  def get_file(self):
    return self.complete_file % (self.plots, self.highlights, self.links)


def get_circos_GC_config_files(biodatabase_name, accession_list):
    from chlamdb.plots import GC
    '''
    accessions: in case of several chromosomes or combinations of chromosomes

    '''

    server, db = manipulate_biosqldb.load_db(biodatabase_name)

    final_gc_var = ''
    final_gc_skew = ''

    for accession in accession_list:
        record = db.lookup(accession=accession)

        #print record
        mean_GC = GC.GC(record.seq)
        #print "mean GC", mean_GC
        circos_gc_var = GC.circos_gc_var(record)
        circos_gc_skew = GC.circos_gc_skew(record)
        final_gc_var += circos_gc_var
        final_gc_skew += circos_gc_skew

        #import pylab
        #import numpy as np
        #cumulated_skew = np.cumsum(values)
        #pylab.plot(cumulated_skew)
        #pylab.plot(values[0:-1])
        #pylab.show()

    return (final_gc_var, final_gc_skew)



if __name__ == '__main__':

    import shell_command
    parser = OptionParser()

    parser.add_option("-i", "--input",dest="input_file",action="store", type="string", help="genbank file", metavar="FILE")
    parser.add_option("-o", "--output",dest="output_file",action="store",default="record", type="string", help="genbank file", metavar="FILE")
    (options, args) = parser.parse_args()



    gb_file = options.input_file
    # get list of all records present in the gbk file
    record_list=[]

    record_list = list(SeqIO.parse(open(gb_file,"r"), "genbank"))
        #print gb_record



    # format output tab delimited table
    #print "contig\ttype\tstart\tstop\tlength\tGC\tstrand\tgene\tfunction\tinference\tgi\tlocus\ttranslation\tsequence"
    #for record in record_list:
    #  for feature in record.features:
    #    if feature is None:
    #      continue
    #    if feature.type == "source":
    #      pass
          #print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t \t " % (feature.contig, feature.type, feature.start, feature.stop, feature.length, feature.GC, feature.strand, feature.gene, feature.product, feature.inference, feature.gi, feature.locus)
    #    else:
    #      print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (feature.contig, feature.type, feature.start, feature.stop, feature.length, feature.GC, feature.strand, feature.gene, feature.product, feature.inference, feature.gi, feature.locus, feature.translation, feature.seq)

    print_circos_record_file(record_list, 'circos_contigs.txt')
    print_circos_gene_file(record_list, strand="-1", out_name = 'circos_minus.txt')
    print_circos_gene_file(record_list, strand="1", out_name = 'circos_plus.txt')
    circos_GC, circos_skew = print_circos_GC_file(record_list)


    circos_conf = Circos_config('circos_contigs.txt')

    circos_conf.add_highlight('circos_plus.txt', fill_color="grey_a1", r1="0.8r", r0="0.75r")
    circos_conf.add_highlight('circos_minus.txt', fill_color="grey_a1", r1="0.75r", r0="0.7r")

    conditions = circos_conf.template_rules % (circos_conf.template_rule('var(value) < 0', 'lred') +
                                                        circos_conf.template_rule('var(value) > 0', 'lblue'))
    circos_conf.add_plot(circos_GC, fill_color="green", r1="0.99r", r0= "0.82r", type="line", rules=conditions)


    t = open('circos_config.txt', "w")
    t.write(circos_conf.get_file())
    t.close()
    cmd = "circos -outputfile circos -outputdir . -conf circos_config.txt"
    #print cmd
    (stdout, stderr, return_code) = shell_command.shell_command(cmd)
