#! /usr/bin/env python

# create tabulated data based on genbank file
# produce one ptt / record present in the genbank file
# TODO: unnecessarily convert BioSeq Records into Record and Feature objects, could be skipped
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------

from Bio import SeqIO
from optparse import OptionParser
import re
from Bio import SeqUtils

def find_index(pattern, seq):
  """Return first item in sequence where f(item) == True."""
  for item in seq:
    if re.match(pattern,item): 
      return seq.index(item)


class Feature:
    def __init__(self):
        
        self.contig = "-" 
        self.type = "-" 
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

class Record:
    def __init__(self, record):
        self.seq = record.seq
        #print "length", len(self.seq)
        try:
          self.contig =  record.name
        except:
          self.contig = "-"
        #print "name", self.contig
        self.features = self.get_one_record_features(record)

    def get_one_record_features(self, one_record):
            
        feature_list = []
        for i in range(0,len(one_record.features)):
            #print one_record.features[i]
            if one_record.features[i].type == "misc_feature":
                continue
            new_feature = Feature() 
            #print  one_record.features[i]
            new_feature.type = one_record.features[i].type
            new_feature.contig = one_record.name
            new_feature.start = one_record.features[i].location.start
            new_feature.stop = one_record.features[i].location.end
            new_feature.length = len(one_record.features[i].location)
            new_feature.strand = one_record.features[i].strand
            try:
                new_feature.gene = one_record.features[i].qualifiers['gene'][0]
            except:
                pass
            try:
                gi_position=find_index("GO*",one_record.features[i].qualifiers['db_xref'])
                new_feature.gi = one_record.features[i].qualifiers['db_xref'][gi_position][3:]
            except:
                pass
         
            #geneID= one_record.features[i].qualifiers['db_xref'][1][7:]
            try:
                new_feature.locus = one_record.features[i].qualifiers['locus_tag'][0]
            except:
                pass
            try:
                new_feature.protein_id = one_record.features[i].qualifiers['protein_id'][0]
            except:
                pass
            try:
                new_feature.product = one_record.features[i].qualifiers['product'][0]
            except:
                pass
            try:
              
              new_feature.inference = one_record.features[i].qualifiers['inference']
            except:
              pass
            try:
                new_feature.translation = one_record.features[i].qualifiers['translation'][0]
            except:
                pass
            new_feature.seq = one_record.features[i].extract(self.seq)
          
            new_feature.GC = SeqUtils.GC(new_feature.seq)
            feature_list.append(new_feature)

        return feature_list


def gbk2table(gb_file_or_record, output_file):
    from Bio.SeqRecord import SeqRecord

    record_list=[]

    if isinstance(gb_file_or_record, list):
        record_list = gb_file_or_record

    elif isinstance(gb_file_or_record, str):
        for gb_record in SeqIO.parse(open(gb_file_or_record,"r"), "genbank") :
            #print gb_record
            record_list.append(Record(gb_record))
    elif isinstance(gb_file_or_record, SeqRecord):
        record_list = [Record(gb_file_or_record)]
    else:
        raise('Unknown input format: provide either path to gbk file or a SeqRecord object of list of SeqRecords')

    with open(output_file, 'w') as f:
        # format output tab delimited table
        f.write("contig\ttype\tstart\tstop\tlength\tGC\tstrand\tgene\tfunction\tinference\tgi\tlocus\ttranslation\tsequence\n")
        for record in record_list:
          for feature in record.features:
            if feature is None:
              continue
            if feature.type == "source":
              pass
              #print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t \t " % (feature.contig, feature.type, feature.start, feature.stop, feature.length, feature.GC, feature.strand, feature.gene, feature.product, feature.inference, feature.gi, feature.locus)
            else:
              f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (feature.contig, feature.type, feature.start, feature.stop, feature.length, feature.GC, feature.strand, feature.gene, feature.product, feature.inference, feature.gi, feature.locus, feature.translation, feature.seq))



if __name__ == '__main__':
    
    
    parser = OptionParser()

    parser.add_option("-i", "--input",dest="input_file",action="store", type="string", help="genbank file", metavar="FILE")
    parser.add_option("-o", "--output",dest="output_file",action="store",default="record", type="string", help="genbank file", metavar="FILE")
    (options, args) = parser.parse_args()

    # get list of all records present in the gbk file
    gbk2table(options.input_file, options.output_file)
