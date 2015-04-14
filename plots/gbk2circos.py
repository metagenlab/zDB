#! /usr/bin/python

# produce one ptt / record present in the genbank file

from Bio import SeqIO
from optparse import OptionParser
import re
from Bio import SeqUtils
import orthogroup_identity_db
import manipulate_biosqldb
import mysqldb_load_mcl_output


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
        return []
    contigs = []
    for i in range(0, len(gap_locations)):
        if i == 0:
            contigs.append([record.name + "_%s" % 1, 0, int(gap_locations[i].start)])
        else:
            contigs.append([record.name + "_%s" % (i+1), int(gap_locations[i-1].end), int(gap_locations[i].start)])
    contigs.append([record.name + "_%s" % (len(gap_locations)+1), int(gap_locations[-1].end), len(record.seq)])
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
              if draft_coordinates:
                  if i%2 == 0:
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], 0, contig[2]-contig[1], 3))
                  else:
                      #print "chr -", record.contig, record.contig, "0",len(record.seq), "spectral-5-div-%s" % (4)
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], 0, contig[2]-contig[1], 4))
                  i+=1
              else:
                  if i%2 == 0:
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], contig[1], contig[2], 3))
                  else:
                      #print "chr -", record.contig, record.contig, "0",len(record.seq), "spectral-5-div-%s" % (4)
                      f.write("chr - %s %s %s %s spectral-5-div-%s\n" % (contig[0], contig[0], contig[1], contig[2], 4))
                  i+=1

          x+=1
      except:
            if i%2 == 0:
                f.write("chr - %s %s 0 %s spectral-5-div-%s\n" % (record.id, record.id, len(record.seq), 3))
            else:
                #print "chr -", record.contig, record.contig, "0",len(record.seq), "spectral-5-div-%s" % (4)
                f.write("chr - %s %s 0 %s spectral-5-div-%s\n" % (record.id, record.id, len(record.seq), 4))
            i += 1
            x += 1
    f.close()

def print_circos_gene_file(record_list, feature_type="CDS", strand ="1",
                           out_name = "circos_genes_plus.txt",
                           locus_highlight=[], draft_data=False,
                           group_id2orthologs_presence=False, query_taxon_id=False,
                           color_missing=True, draft_coordinates=False):

    print "highlight:", locus_highlight
    
    import numpy
    print "draft_data", draft_data
    if strand == "1":
        f = open(out_name, "w")
    if strand == "-1":
        f = open(out_name, "w")


    for y, record in enumerate(record_list):
        for feature in record.features:
            if feature.type == feature_type:
                if str(feature.strand) == strand:
                    try:
                        #print "draft_data[y]", draft_data[y]
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
                    # in case the second record is not fragmented in several contigs (no draft data)
                    except IndexError:
                        print "no draft for", record.id
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    except TypeError:
                        print "no draft for", record.id
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end
                    #print "longueur", len(feature.location), type(len(feature.location)), feature.location.start, feature.location.end

                    if numpy.abs(feature.location.start-feature.location.end) > 50000:
                        continue

                    if 'pseudo' in feature.qualifiers:

                        continue
                    try:
                        a = feature.qualifiers['orthogroup']
                    except:
                        print "problem with", feature
                        continue

                    if feature.qualifiers['orthogroup'][0] in locus_highlight:

                        print "COLORS!!!!!!!!!!!!!!!!!!!"
                        try:
                            '''
                            f.write('%s %s %s fill_color=violet, id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                                     start,
                                                                                                     end,
                                                                                                     feature.qualifiers['locus_tag'][0],
                                                                                                     feature.qualifiers['gene'][0],
                                                                                                     re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                                     feature.qualifiers['product'][0])))
                            '''
                            f.write('%s %s %s fill_color=violet, id=%s\n' % (contig,
                                                                    start,
                                                                    end, feature.qualifiers['locus_tag'][0]))


                        except:
                            '''
                            f.write('%s %s %s fill_color=violet, id=locus:_%s_gene:_%s_product:_%s\n' % (contig,
                                                                                                     start,
                                                                                                     end,
                                                                                                     feature.qualifiers['locus_tag'][0],
                                                                                                     "-",
                                                                                                     re.sub("[ |\-|(|)|\[|\]|\.|,]", "_",
                                                                                                     feature.qualifiers['product'][0])))
                            '''
                            f.write('%s %s %s fill_color=violet, id=%s\n' % (contig,
                                                                    start,
                                                                    end,
                                                                    feature.qualifiers['locus_tag'][0]))



                    elif group_id2orthologs_presence and query_taxon_id and color_missing:
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
                if str(feature.strand) == strand:
                    try:
                        for i in draft_data[y]:
                            if draft_coordinates:
                                if feature.location.start >= i[1] and feature.location.end <= i[2]:
                                    contig = i[0]
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    contig = i[0]
                                    start = feature.location.start
                                    end = feature.location.end
                    except:
                        contig = record.id # fill_color=violet
                        start = feature.location.start
                        end = feature.location.end


                    f.write('%s %s %s fill_color=pblue\n ' % (contig,
                                                            start,
                                                            end))


            if feature.type == 'tRNA':
                if str(feature.strand) == strand:
                    try:
                        for i in draft_data[y]:
                            if feature.location.start >= i[1] and feature.location.end <= i[2]:
                                contig = i[0]
                                start = feature.location.start - i[1]
                                end = feature.location.end - i[1]
                    except:
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




def print_circos_GC_file(record_list, feature_type="CDS", out_name="circos_GC.txt", draft_data = False):


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
    import GC

    f = open('circos_GC_var_%s.txt' % record_list[0].name, 'w')
    g = open('circos_GC_skew_%s.txt' % record_list[0].name, 'w')

    out_var = ''
    out_skew = ''
    for record in record_list:
        out_var += GC.circos_gc_var(record)
        out_skew += GC.circos_gc_skew(record)
    f.write(out_var)
    g.write(out_skew)





  
def orthology_circos_files(server, record_list, reference_taxon_id, biodatabase_name, out_dir, locus_highlight=[],
                           feature_type="CDS", taxon_list=False, draft_data=False, query_taxon_id=False,
                           draft_coordinates=False, color_missing=True):

    print "highlight", locus_highlight

    import os
    print "orthology_circos_files"
    print "draft_data", draft_data
    print "locus highlight", locus_highlight
    all_file_names = {}
    ortho_size = mysqldb_load_mcl_output.get_all_orthogroup_size(server, biodatabase_name)
    ortho_identity = orthogroup_identity_db.orthogroup2average_identity(biodatabase_name)
    # get dictionnary orthogroup_id2family_size
    ortho_family_size = mysqldb_load_mcl_output.get_family_size(server, biodatabase_name)


    print ortho_family_size

    ortho_table = "orthology_%s" % biodatabase_name
    print "ortho_table", ortho_table
    group_id2orthologs_presence = mysqldb_load_mcl_output.get_orthology_matrix_merging_plasmids(server, biodatabase_name)
    # get taxons ids used as column names in orthology table
    if taxon_list is False:
        taxon_list = manipulate_biosqldb.get_column_names(server, ortho_table)[1:]

    #print "all_taxon_ids", taxon_list

    taxon_id2destription = {}

    for taxon in taxon_list:
        sql = 'select accession from bioentry where taxon_id=%s limit 1;' % taxon
        description = server.adaptor.execute_and_fetchall(sql)[0]
        taxon_id2destription[taxon] = description[0]



    #print "taxon_id2destription", taxon_id2destription


    all_file_names["contigs"] = os.path.join(out_dir,
                                           "circos_contigs_%s.txt" % taxon_id2destription[reference_taxon_id])
    all_file_names["minus"] = os.path.join(out_dir,
                                         "circos_genes_minus_%s.txt" % taxon_id2destription[reference_taxon_id])
    all_file_names["plus"] = os.path.join(out_dir,
                                        "circos_genes_plus_%s.txt" % taxon_id2destription[reference_taxon_id])
    all_file_names["GC"] = os.path.join(out_dir,
                                      "circos_GC_%s.txt" % taxon_id2destription[reference_taxon_id])
    all_file_names["orthogroups"] = os.path.join(out_dir,
                                               "circos_orthogroup_size_%s.txt"  % taxon_id2destription[reference_taxon_id])
    all_file_names["orthogroups_identity"] = os.path.join(out_dir,
                                               "circos_orthogroup_identity_%s.txt"  % taxon_id2destription[reference_taxon_id])
    all_file_names["orthogroups_family"] = os.path.join(out_dir,
                                               "circos_orthogroup_family_size_%s.txt"  % taxon_id2destription[reference_taxon_id])

    print "writing record file"
    print_circos_record_file(record_list, out_name = all_file_names["contigs"], draft_data=draft_data)
    print "writing minus strand file"
    print_circos_gene_file(record_list, strand="-1", out_name = all_file_names["minus"], draft_data=draft_data,
                         locus_highlight=locus_highlight, group_id2orthologs_presence=group_id2orthologs_presence,
                         query_taxon_id=query_taxon_id)
    print "writing plus strand file"
    print_circos_gene_file(record_list, strand="1", out_name = all_file_names["plus"], draft_data=draft_data,
                         locus_highlight=locus_highlight, group_id2orthologs_presence=group_id2orthologs_presence,
                         query_taxon_id=query_taxon_id)
    print "writing GC file"
    print_circos_GC_file(record_list, out_name = all_file_names["GC"], draft_data=draft_data)




    f = open(all_file_names["orthogroups"], "w")

    g = open(all_file_names["orthogroups_identity"], "w")
    h = open(all_file_names["orthogroups_family"], "w")
    taxon_files = []
    all_file_names["genomes"] = []

    #print "all_taxon_ids2", taxon_list


    for i in range(0, len(taxon_list)):
        print i
        #print "taxon", taxon_list[i]
        if int(taxon_list[i]) == int(reference_taxon_id):
            rem = i
        else:
            file_name = os.path.join(out_dir, "circos_taxon_%s_vs_%s.txt" % (taxon_id2destription[reference_taxon_id], taxon_id2destription[taxon_list[i]]))
            all_file_names["genomes"].append(file_name)
            taxon_files.append(open(file_name, "w"))

    taxon_list.pop(rem)
    y = 0
    for record in record_list:
        for feature in record.features:
            import time
            #time.sleep(2)
            if feature.type == feature_type:
                if not 'pseudo' in feature.qualifiers:
                    try:
                        for i in draft_data[y]:
                            #print "TZDHGFHFDHD,", i, int(feature.location.start), int(feature.location.end), record.id
                            if int(feature.location.start) >= i[1] and int(feature.location.end) <= i[2]:
                                #print "OKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK"
                                contig = i[0]
                                if draft_coordinates:
                                    start = feature.location.start - i[1]
                                    end = feature.location.end - i[1]
                                else:
                                    start = feature.location.start
                                    end = feature.location.end

                    except:
                        contig = record.id
                        start = feature.location.start
                        end = feature.location.end
                    #print "contig", contig
                    try:
                        line = "%s %s %s %s\n" % (contig, start, end, ortho_size[feature.qualifiers['orthogroup'][0]])
                    except:
                        print "problem ith", feature
                        continue
                    line2 = "%s %s %s %s\n" % (contig, start, end, ortho_identity[feature.qualifiers['orthogroup'][0]])
                    line3 = "%s %s %s %s\n" % (contig, start, end, ortho_family_size[feature.qualifiers['orthogroup'][0]])
                    f.write(line)
                    g.write(line2)
                    h.write(line3)
                    for taxon_id, taxon_file in zip(taxon_list, taxon_files):
                        #print group_id2orthologs_presence[feature.qualifiers['orthogroup'][0]]
                        if group_id2orthologs_presence[feature.qualifiers['orthogroup'][0]][int(taxon_id)] > 0:
                            if feature.qualifiers['orthogroup'][0] not in locus_highlight:
                                taxon_file.write("%s %s %s id=%s\n" % (contig, start, end, feature.qualifiers['locus_tag'][0]))
                            else:
                                print "COLORS!!!!!!!!!!!!!!!!!!!"
                                taxon_file.write("%s %s %s fill_color=violet, id=%s\n" % (contig, start, end, feature.qualifiers['locus_tag'][0]))
                        elif color_missing:
                            taxon_file.write("%s %s %s fill_color=black\n" % (contig, start, end))
                        else:
                            pass #taxon_file.write("%s %s %s\n" % (contig, feature.start, feature.stop))
        y += 1
    return (all_file_names, taxon_id2destription)
        
  #print manipulate_biosqldb.get_taxon_id_list(server, biodatabase_name)


class Circos_config:
  def __init__(self, caryotype_file, chr_spacing_list=[], show_ticks="yes", show_tick_labels="yes", ideogram_spacing=0, radius=0.75, label_radius=0.175):
    import re

    self.plots = ""
    self.links = ""
    self.highlights = ""
    
    self.template_caryotype= "karyotype = %s\n" \
                             " chromosomes_units           = 1000\n" \
                             " chromosomes_display_default = yes\n" % caryotype_file


    self.template_ideograms = "<ideogram>\n" \
                              " <spacing>\n" \
                              " default            = %su\n" \
                              " %s" \
                              " </spacing>\n" \
                              " \n" \
                              " # thickness and color of ideograms\n" \
                              " thickness          = 15p\n" \
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
                              " show_label         = yes\n" \
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
                              " </ideogram>\n" % (ideogram_spacing, self.add_spacing(chr_spacing_list), radius, label_radius)

    self.template_ticks = "show_ticks         = %s\n" \
                          " show_tick_labels   = %s\n" \
                          " \n" \
                          " <ticks>\n" \
                          " tick_label_font    = condensed\n" \
                          " radius             = dims(ideogram,radius_outer)\n" \
                          " label_offset       = 8p\n" \
                          " label_size         = 8p\n" \
                          " color              = black\n" \
                          " thickness          = 4p\n" \
                          " \n" \
                          " <tick>\n" \
                          " spacing           = 100u\n" \
                          " size              = 16p\n" \
                          " label_multiplier  = 1e-3\n" \
                          " show_label        = yes\n" \
                          " label_size        = 35p\n" \
                          " format            = %s kb\n" \
                          " thickness         = 5p\n" \
                          " </tick>\n" \
                          " \n" \
                          " <tick>\n" \
                          " spacing           = 10u\n" \
                          " size              = 8p\n" \
                          " show_label        = no\n" \
                          " label_size        = 5p\n" \
                          " format            = %s\n" \
                          " </tick>\n" \
                          " \n" \
                          " </ticks>\n" % (show_ticks, show_tick_labels, "%%d", "%%d")
    print self.template_ticks




    self.template_rules = "<rules>\n" \
                          " %s\n" \
                          "</rules>\n"


    self.settings ="<colors>\n" \
                   " <<include colors.rn.conf>>\n" \
                   " <<include brewer.all.conf>>\n" \
                   " </colors>\n" \
                   " <image>\n"\
                   " image_map_use      = yes\n" \
                   " image_map_overlay  = yes\n" \
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

    self.complete_file = self.template_caryotype + self.template_ideograms + self.template_ticks + "%s %s %s" + self.settings

  def _template_spacing(self, chr1, chr2):
    template = '<pairwise %s %s>\n' \
               ' spacing = 8u\n' \
               '</pairwise>\n' % (chr1, chr2)
    return template
    
  def _template_plot(self, file, type="line", r0=1, r1=1.05, color="black", fill_color="red", thickness = "2p", z = 1, rules =""):
    template = "<plot>\n" \
               "type		    = %s\n" \
               " #min               = -0.4\n" \
               " #max               = 0.4\n" \
               " r0                 = %s\n" \
               " r1                 = %s\n" \
               " color              = %s\n" \
               " fill_color         = %s\n" \
               " thickness          = %s\n" \
               " file               = %s\n" \
               " z                  = %s\n" \
               " %s\n" \
               " </plot>\n" % (type, r0, r1, color, fill_color, thickness, file, z, rules)
    return template

  def template_rule(self, condition, fill_color):
    template = "<rule>\n" \
               " condition          = %s\n" \
               " fill_color         = %s\n" \
               " </rule>\n" % (condition, fill_color)
    return template

  def _template_link(self, link_file, color="black_a5", thickness=1):
    template = "<link>\n" \
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


  def add_plot(self, file, type="line", r0=1, r1=1.05, color="black", fill_color="grey_a1", thickness = "2p", z = 1, rules =""):
    plot = self._template_plot(file, type, r0, r1, color, fill_color, thickness, z, rules)
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


  def add_highlight(self, file, fill_color="grey_a1", r1="1.55r", r0="1.50r"):
    highlight = self._template_highlight(file, fill_color, r1, r0)
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
    print self.complete_file
    return self.complete_file % (self.plots, self.highlights, self.links)


def get_circos_GC_config_files(biodatabase_name, accession_list):
    import GC
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
        print "mean GC", mean_GC
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
    
    
    parser = OptionParser()

    parser.add_option("-i", "--input",dest="input_file",action="store", type="string", help="genbank file", metavar="FILE")
    parser.add_option("-o", "--output",dest="output_file",action="store",default="record", type="string", help="genbank file", metavar="FILE")
    (options, args) = parser.parse_args()



    gb_file = options.input_file
    # get list of all records present in the gbk file
    record_list=[]
   
    for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
        #print gb_record
        record_list.append(Record(gb_record))
    

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
        
    print_circos_record_file(record_list)
    print_circos_gene_file(record_list)
    print_circos_GC_file(record_list)
    print_circos_GC_file(record_list)
