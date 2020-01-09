from Bio import Entrez, SeqIO
from Bio.SeqUtils import CheckSum
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

import pandas as pd
import sys
import sqlite3
import urllib
import gzip
import re
import os
import datetime
import logging
import http.client

from collections import Iterable

def merge_gbk(gbk_records, filter_size=0, gi=False):
    '''
    merge multiple contigs into a single DNA molecule with 200*N between contigs
    keep source description from the first record
    remove contigs smaller than <filter_size>

    :param gbk_records:
    :param filter_size:
    :param gi:
    :return:
    '''

    n=0
    if len(gbk_records) == 1:
        merged_rec = gbk_records[0]
    else:
        for i, rec in enumerate(gbk_records):
            # remove source feature of all records except the first one
            if rec.features[0].type == 'source' and i != 0:
                rec.features.pop(0)
            # filter small contigs
            if len(rec) > filter_size:
                if n == 0:
                    n+=1
                    merged_rec = rec
                else:
                    merged_rec+=rec
                # you could insert a spacer if needed
                # do not add spacer after the last contig
                if i != len(gbk_records)-1:
                    merged_rec += "N" * 200

                    my_start_pos = ExactPosition(len(merged_rec)-200)
                    my_end_pos = ExactPosition(len(merged_rec))
                    my_feature_location = FeatureLocation(my_start_pos, my_end_pos)
                    my_feature = SeqFeature(my_feature_location, type="assembly_gap")
                    merged_rec.features.append(my_feature)

    try:
        merged_rec.id = gbk_records[0].annotations["accessions"][-1]
    except KeyError:
        merged_rec.id = gbk_records[0].id

    if gi:
        merged_rec.annotations["gi"] = gi

    merged_rec.description = "%s" % gbk_records[0].annotations["organism"]
    merged_rec.annotations = gbk_records[0].annotations                                             
    try:
        merged_rec.name = gbk_records[0].annotations["accessions"][-1]
    except KeyError:
        merged_rec.name = gbk_records[0].id
    my_start_pos = ExactPosition(0)
    my_end_pos = ExactPosition(len(merged_rec))
    merged_rec.features[0].location = FeatureLocation(my_start_pos, my_end_pos)

    return merged_rec



def filter_plasmid(record_list):
    plasmid_record_list = []
    chromosome_record_list = []

    for record in record_list:
        # plasmid.annotations['organism']
        if record.features[0].type == 'source':
            if 'plasmid' in record.description or "plasmid" in record.features[0].qualifiers:
                plasmid_record_list.append(record)
            else:
                chromosome_record_list.append(record)
        else:
            if 'plasmid' in record.description:
                plasmid_record_list.append(record)
            else:
                chromosome_record_list.append(record)
    return (chromosome_record_list, plasmid_record_list)


def count_missing_locus_tags(gbk_record):
    count_CDS = 0
    count_no_locus = 0
    for feature in gbk_record.features:
        if feature.type == 'CDS':
            count_CDS += 1
            if "locus_tag" not in feature.qualifiers:
                count_no_locus += 1
    return count_no_locus, count_CDS

def is_annotated(gbk_record):
    return not (len(gbk_record.features) == 1
            and gbk_record.features[0].type == 'source')

def update_record_taxon_id(record, n):
    if record.features[0].type == 'source':
        if 'db_xref' in record.features[0].qualifiers:
            for item in record.features[0].qualifiers['db_xref']:
                if 'taxon' in item:
                    index = record.features[0].qualifiers['db_xref'].index(item)
                    record.features[0].qualifiers['db_xref'][index] = "taxon:%s" % n
    else:
        print('ACHRTUNG\t no source for record \t%s' % record.name)
    return record


def rename_source(record):
    if 'strain' in record.features[0].qualifiers:

        print('--', record.features[0].qualifiers['strain'][0])
        if ';' in record.features[0].qualifiers['strain'][0]:
            print('ACHRTUNG: record has 2 strain names! \t%s\t --> check and edit source manually' % record.name)
            # put everythink lower size
            strain = record.features[0].qualifiers['strain'][0].split(';')[1]
        else:
            strain = record.features[0].qualifiers['strain'][0]
        if strain == 'strain':
            return (False, False)
        if strain.lower() not in record.annotations['source'].lower():
            msg = '%s' % record.annotations['source']
            print("ACHTUNG changing source\t%s\t--> %s " % (msg, record.annotations['source'] + strain))


        return strain, "%s %s" % (record.annotations['source'], strain)
    else:
        return (False, False)

def clean_description(description):
    import re
    description = re.sub(", complete genome\.", "", description)
    description = re.sub(", complete genome", "", description)
    description = re.sub(", complete sequence\.", "", description)
    description = re.sub("strain ", "", description)
    description = re.sub("str\. ", "", description)
    description = re.sub(" complete genome sequence\.", "", description)
    description = re.sub(" complete genome\.", "", description)
    description = re.sub(" chromosome", "", description)
    description = re.sub(" DNA", "", description)
    description = re.sub("Merged record from ", "", description)
    description = re.sub(", wgs", "", description)
    description = re.sub("Candidatus ", "", description)
    description = re.sub(".contig.0_1, whole genome shotgun sequence.", "", description)
    description = re.sub("complete genome, isolate", "", description)
    description = re.sub(" complete", "", description)
    description = re.sub(" genome assembly.*", "", description)
    description = re.sub("Chlamydophila", "Chlamydia", description)
    description = re.sub(", whole genome shotgun sequence", "", description)

    return description


def check_gbk(gbff_files, minimal_contig_length=1000):
    reannotation_list = []

    # count the number of identical source names
    source2count = {}
    accession2count = {}
    for genome_number, gbff_file in enumerate(gbff_files):

        records = list(SeqIO.parse(gzip.open(gbff_file, "rt"), "genbank"))

        for record in records:
            n_missing, total = count_missing_locus_tags(record)
            if n_missing > 0:
                print ('Warrning: %s/%s missing locus tag for record %s' % (n_missing, total, record.name))

        chromosome, plasmids = filter_plasmid(records)

        cleaned_records = []
        plasmid_reannot = False
        chromosome_reannot = False

        for n_plasmid, plasmid in plasmids:
            annot_check = is_annotated(plasmid)
            if annot_check:

                plasmid.description = clean_description(plasmid.description)

                plasmid = update_record_taxon_id(plasmid, 1000 + genome_number)
                strain, new_source = rename_source(plasmid)
                print("plasmid:", strain, new_source )
                if new_source:
                    if not 'plasmid' in new_source:
                        new_source = "%s plasmid %s" % (new_source, n_plasmid+1)
                    if strain.lower() not in plasmid.annotations['source'].lower():
                        plasmid.description = new_source
                    if strain.lower() not in plasmid.annotations['organism'].lower():
                        plasmid.annotations['organism'] = new_source
                    if strain.lower() not in plasmid.annotations['source'].lower():
                        plasmid.annotations['source'] = new_source
                else:
                    print ('ACHTUNG\t no strain name for \t%s\t, SOUCE uniqueness should be checked manually' % merged_record.id)
                # check if accession is meaningful
                if 'NODE_' in plasmid.id or 'NODE_' in plasmid.name:
                    print ('ACHTUNG\t accession probably not unique (%s) for \t%s\t --> should be checked manually' % (merged_record.id))
                cleaned_records.append(plasmid)
            else:
                plasmid_reannot = True
                print("Warrning: unannotated genome: %s" % plasmid)

        if len(chromosome) > 0:

            '''
            Assume single chromosome bacteria.
            If multiple record founds, consider those as contigs.
            Contigs contatenation with 200 N between each (labelled as assembly_gap feature)
            '''

            ########## chromosome ###########
            if chromosome[0].seq == 'N'*len(chromosome[0].seq):
                print('Warning: No sequences for %s, skipping' % gbff_file)
                continue
            annot_check = is_annotated(chromosome[0])
            if annot_check:
                if len(chromosome) > 1:
                    merged_record = merge_gbk(chromosome)
                else:
                    merged_record = chromosome[0]
                print(merged_record.description )
                merged_record.description = clean_description(merged_record.description)
                print(merged_record.description)
                # rename source with strain name
                merged_record = update_record_taxon_id(merged_record, 1000 + genome_number)
                strain, new_source = rename_source(merged_record)
                if new_source:
                    if strain.lower() not in merged_record.annotations['source'].lower():
                        merged_record.description = new_source
                    if strain.lower() not in merged_record.annotations['organism'].lower():
                        merged_record.annotations['organism'] = new_source
                    if strain.lower() not in merged_record.annotations['source'].lower():
                        merged_record.annotations['source'] = new_source
                else:
                    print('ACHTUNG\t no strain name for\t%s' % gbff_file)
                # check if accession is meaningful
                if 'NODE_' in merged_record.id or 'NODE_' in merged_record.name:
                    print('ACHTUNG\t accession probably not unique (%s) for\t%s' % (merged_record.id, gbff_file))
                cleaned_records.append(merged_record)
            else:
                chromosome_reannot = True
                print("Warrning: unannotated genome: %s" % chromosome)


        if plasmid_reannot and not chromosome_reannot and len(chromosome) > 0:
            print("plasmid", plasmid_reannot)
            print("chr", chromosome_reannot)
            raise TypeError('Combination of unannotated plasmid(s) and annotated chromosome')
        elif not plasmid_reannot and chromosome_reannot and len(plasmids) > 0:
            print("plasmid", plasmid_reannot)
            print("chr", chromosome_reannot)
            raise TypeError('Combination of annotated plasmid(s) and unannotated chromosome')
        elif plasmid_reannot or chromosome_reannot:
            print("plasmid", plasmid_reannot)
            print("chr", chromosome_reannot)
            raise TypeError('Some genomes are not annotated!')
        else:
            out_name = gbff_file.split('.')[0] + '_merged.gbk'
            with open(out_name, 'w') as f:
                SeqIO.write(cleaned_records, f, 'genbank')

def filter_out_unannotated(gbk_file):
    # assumes ".gbk" file ending
    result_file = gbk_file[:-len(".gbk")] + "_filtered.gbk"
    to_keep = []
    for record in SeqIO.parse(gbk_file, "genbank"):
        if is_annotated(record):
            to_keep.append(record)

    SeqIO.write(to_keep, result_file, "genbank")

def string_id2pubmed_id_list(accession):
    link = 'http://string-db.org/api/tsv/abstractsList?identifiers=%s' % accession
    try:
        data = urllib2.urlopen(link).read().rstrip().decode('utf-8').split('\\n')[1:]
    except urllib2.URLError:
        return False
    pid_list = [row.split(':')[1] for row in data]
    return pid_list

def get_string_PMID_mapping(string_map):
    o = open("string_mapping_PMID.tab", "w")
    with open(string_map, 'r') as f:
        for n, row in enumerate(f):
            if n == 0:
                continue
            data = row.rstrip().split("\t")
            pmid_list = string_id2pubmed_id_list(data[1])
            if pmid_list:
                for id in pmid_list:
                    o.write("%s\t%s\n" % (data[0], id))
            else:
                o.write("%s\tNone\n" % (data[0]))


def get_pdb_mapping(fasta_file, database_dir):
    conn = sqlite3.connect(database_dir + "/pdb/pdb.db")
    cursor = conn.cursor()
    pdb_map = open('pdb_mapping.tab', 'w')
    no_pdb_mapping = open('no_pdb_mapping.faa', 'w')

    pdb_map.write("locus_tag\tpdb_id\n")

    records = SeqIO.parse(fasta_file, "fasta")
    no_pdb_mapping_records = []
    for record in records:
        sql = 'select accession from hash_table where sequence_hash=?'
        cursor.execute(sql, (CheckSum.seguid(record.seq),))
        hits = cursor.fetchall()
        if len(hits) == 0:
            no_pdb_mapping_records.append(record)
        else:
            for hit in hits:
              pdb_map.write("%s\t%s\n" % (record.id,
                                              hit[0]))
    SeqIO.write(no_pdb_mapping_records, no_pdb_mapping, "fasta")


def get_tcdb_mapping(fasta_file, database_dir):
    conn = sqlite3.connect(database_dir + "/TCDB/tcdb.db")
    cursor = conn.cursor()
    tcdb_map = open('tcdb_mapping.tab', 'w')
    no_tcdb_mapping = open('no_tcdb_mapping.faa', 'w') 
    tcdb_map.write("locus_tag\ttcdb_id\n")

    records = SeqIO.parse(fasta_file, "fasta")
    no_tcdb_mapping_records = []
    for record in records:
        sql = 'select accession from hash_table where sequence_hash=?'
        cursor.execute(sql, (CheckSum.seguid(record.seq),))
        hits = cursor.fetchall()
        if len(hits) == 0:
            no_tcdb_mapping_records.append(record)
        else:
            for hit in hits:
              tcdb_map.write("%s\t%s\n" % (record.id,
                                              hit[0]))

    SeqIO.write(no_tcdb_mapping_records, no_tcdb_mapping, "fasta")


def get_string_mapping(fasta_file, database_dir):
    conn = sqlite3.connect(database_dir + "/string/string_proteins.db")
    cursor = conn.cursor()
    string_map = open('string_mapping.tab', 'w')
    string_map.write("locus_tag\tstring_id\n")

    records = SeqIO.parse(fasta_file, "fasta")
    no_mapping_string_records = []
    for record in records:
        sql = 'select accession from hash_table where sequence_hash=?'
        cursor.execute(sql, (CheckSum.seguid(record.seq),))
        hits = cursor.fetchall()
        for hit in hits:
          string_map.write("%s\t%s\n" % (record.id, hit[0]))


# Maybe possible to compact this one with 
# get_nr_sequences to have less lines of code and 
# be more efficient
# TODO : check with Trestan
def convert_gbk_to_faa(gbf_file, edited_gbf):
    records = SeqIO.parse(gbf_file, 'genbank')
    edited_records = open(edited_gbf, 'w')

    for record in records:
        protein2count = {}
        for feature in record.features:
            if (feature.type == 'CDS'
                    and 'pseudo' not in feature.qualifiers
                    and 'pseudogene' not in feature.qualifiers):

                if "locus_tag" in feature.qualifiers:
                    locus_tag = feature.qualifiers["locus_tag"][0]
                else:
                    protein_id = feature.qualifiers["protein_id"][0].split(".")[0]
                    if protein_id not in protein2count:
                        protein2count[protein_id] = 1
                        locus_tag = protein_id
                    else:
                        protein2count[protein_id] += 1
                        locus_tag = "%s_%s" % (protein_id, protein2count[protein_id])
                try:
                    edited_records.write(">%s %s\n%s\n" % (locus_tag,
                                                     record.description,
                                                     feature.qualifiers['translation'][0]))
                except KeyError:
                    print("problem with feature:", feature)

# filter out small sequences and ambiguous amino-acids
def filter_sequences(fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    processed_records = []
    for record in records:
        if len(record.seq) >= 10:
            processed_records.append(SeqRecord(Seq(re.sub("B|Z|J", "X", str(record.seq)), IUPAC.protein),
                                  id=record.id, 
                                  name=record.name,
                                  description=record.description))

    SeqIO.write(processed_records, "filtered_sequences.faa", "fasta")

def get_nr_sequences(fasta_file, genomes_list):
    locus2genome = {}
    for fasta in genomes_list:
        genome = os.path.basename(fasta).split('.')[0]
        for seq in SeqIO.parse(fasta, "fasta"):
            locus2genome[seq.name] = genome
    nr_fasta = open('nr.faa', 'w')
    nr_mapping = open('nr_mapping.tab', 'w')

    checksum_nr_list = []

    records = SeqIO.parse(fasta_file, "fasta")
    updated_records = []

    for record in records:

        checksum = CheckSum.crc64(record.seq)
        nr_mapping.write("%s\t%s\t%s\n" % (record.id,
                                          checksum,
                                          locus2genome[record.id]))
        if checksum not in checksum_nr_list:
            checksum_nr_list.append(checksum)
            record.id = checksum
            record.name = ""
            updated_records.append(record)

    SeqIO.write(updated_records, nr_fasta, "fasta")

def setup_orthology_db(fasta_file, nr_mapping_file, orthogroup):
  fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

  conn = sqlite3.connect("orthology.db")
  cursor = conn.cursor()

  # sequence table
  sql0 = 'create table sequence_hash2aa_sequence (sequence_hash binary, sequence TEXT )'
  cursor.execute(sql0,)
  sql = 'insert into  sequence_hash2aa_sequence values (?, ?)'
  for hash in fasta_dict:
    cursor.execute(sql, (hash, str(fasta_dict[hash].seq)))

  # hash mapping table
  sql1 = 'create table locus_tag2sequence_hash (locus_tag varchar(200), sequence_hash binary)'
  cursor.execute(sql1,)

  sql = 'insert into locus_tag2sequence_hash values (?,?)'
  with open(nr_mapping_file, 'r') as f:
      for row in f:
          data = row.rstrip().split("\t")[0:2]
          cursor.execute(sql, data)
  conn.commit()

  # orthogroup table
  sql2 = 'create table locus_tag2orthogroup (locus_tag varchar(200), orthogroup varchar(200))'
  cursor.execute(sql2,)
  sql = 'insert into locus_tag2orthogroup values (?, ?)'
  with open(orthogroup, 'r') as f:
      for row in f:
          data = row.rstrip().split(" ")
          for locus in data[1:]:
            cursor.execute(sql,(locus, data[0][0:-1]))
  conn.commit()

  # index hash, locus and orthogroup columns
  sql_index_1 = 'create index hash1 on sequence_hash2aa_sequence (sequence_hash);'
  sql_index_2 = 'create index hash2 on locus_tag2sequence_hash (sequence_hash);'
  sql_index_3 = 'create index locus1 on locus_tag2sequence_hash (locus_tag);'
  sql_index_4 = 'create index locus2 on locus_tag2orthogroup (locus_tag);'
  sql_index_5 = 'create index og on locus_tag2orthogroup (orthogroup);'

  cursor.execute(sql_index_1)
  cursor.execute(sql_index_2)
  cursor.execute(sql_index_3)
  cursor.execute(sql_index_4)
  cursor.execute(sql_index_5)
  conn.commit()


# This function parses the result of orthofinder/orthoMCL, 
# retrieves the groups of orthologs detected by the tools
# and returns the core orthologs (i.e. the ortholog present
# in all samples)
#
# Orthofinder output file : Orthogroup_ID: locus1 locus2... locusN
def orthofinder2core_groups(fasta_list,
                          mcl_file,
                          n_missing=0,
                          orthomcl=False):
    orthogroup2locus_list = {}
    with open(mcl_file, 'r') as f:
        for line in f:
            if orthomcl:
                # not sure this code will work with orthoMCL, to be checked
                groups = line.rstrip().split('\t')
                groups = [i.split('|')[1] for i in groups]
            else:
                groups = line.rstrip().split(' ')
                group_id = groups[0][:-1] # remove lagging ':'
                groups = groups[1:]
            orthogroup2locus_list[group_id] = groups

    locus2genome = {}
    for fasta in fasta_list:
        genome = os.path.basename(fasta).split('.')[0]
        for seq in SeqIO.parse(fasta, "fasta"):
            locus2genome[seq.name] = genome

    df = pd.DataFrame(index=orthogroup2locus_list.keys(), columns=set(locus2genome.values()))
    df = df.fillna(0)

    for group_id,loci_list in orthogroup2locus_list.items():
        for locus in loci_list:
            genome = locus2genome[locus]
            df.loc[group_id, genome] += 1

    # why is this necessary? entries are already 
    # numeric values
    df =df.apply(pd.to_numeric, args=('coerce',))

    n_genomes = len(set(locus2genome.values()))
    n_minimum_genomes = n_genomes-n_missing
    freq_missing = (n_genomes-float(n_missing))/n_genomes
    limit = freq_missing*n_genomes

    # Paralogs will have several copies per genomes (df > 1)
    # Why are they removed?
    groups_with_paralogs = df[(df > 1).sum(axis=1) > 0].index
    df = df.drop(groups_with_paralogs)

    # should be an equality?
    core_groups = df[(df == 1).sum(axis=1) >= limit].index.tolist()

    return core_groups, orthogroup2locus_list, locus2genome

def get_core_orthogroups(genomes_list, int_core_missing):
    core_groups, orthogroup2locus_list, locus2genome = orthofinder2core_groups(genomes_list,
          'Orthogroups.txt', int_core_missing, False)

    for group_id in core_groups:
    # sequence_data = SeqIO.to_dict(SeqIO.parse("OG{0:07d}_mafft.faa".format(int(one_group.split('_')[1])), "fasta"))
        sequence_data = SeqIO.to_dict(SeqIO.parse(group_id + "_mafft.faa", "fasta"))
        dest = group_id + '_taxon_ids.faa'
        new_fasta = []
        for locus in orthogroup2locus_list[group_id]:
            tmp_seq = sequence_data[locus]
            tmp_seq.name = locus2genome[locus]
            tmp_seq.id = locus2genome[locus]
            tmp_seq.description = locus2genome[locus]
            new_fasta.append(tmp_seq)

        out_handle = open(dest, 'w')
        SeqIO.write(new_fasta, out_handle, "fasta")
        out_handle.close()

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]

def accession2taxid_entrez(accession):
    from Bio import Entrez
    from socket import error as SocketError
    import errno

    Entrez.email = "trestan.pillonel@chuv.ch"
    Entrez.api_key = "719f6e482d4cdfa315f8d525843c02659408"
    try:
        handle = Entrez.esummary(db="protein", id="%s" % accession, retmax=1)
    except SocketError as e:
        if e.errno != errno.ECONNRESET:
            raise('error connexion with %s' % accession)
        else:
            import time
            print ('connexion error, trying again...')
            time.sleep(60)
            accession2taxid_entrez(accession)
    record = Entrez.read(handle, validate=False)[0]
    # record['AccessionVersion'],
    # record['TaxId'],
    # record['Title'],
    # record['Length']
    return int(record['TaxId'])

def get_refseq_hits_taxonomy(hit_table, database_dir):
    conn = sqlite3.connect("refseq_taxonomy.db")
    cursor = conn.cursor()
    cursor.execute("PRAGMA synchronous = OFF")
    cursor.execute("BEGIN TRANSACTION")
    conn_refseq = sqlite3.connect(database_dir + "/refseq/merged_refseq.db")
    cursor_refseq = conn_refseq.cursor()

    conn_taxid = sqlite3.connect(database_dir + "/ncbi-taxonomy/prot_accession2taxid.db")
    cursor_taxid = conn_taxid.cursor() 
    sql = 'create table refseq_hits (accession varchar(200), taxid INTEGER, description TEXT, length INTEGER)'
    cursor.execute(sql,)
    conn.commit()

    f = open(hit_table, 'r')
    # get list of accessions, remove version number
    hit_list = [i.rstrip().split(".")[0] for i in f]

    chunk_lists = chunks(hit_list, 5000)
    sql_template = 'insert into refseq_hits values (?,?,?,?)'
    hits = []

    # TODO : make a single query using join on the two tables would probably
    # be more efficient and spare some lines of code
    template_annotation = 'select accession, description, sequence_length from refseq where accession in ("%s")'
    template_taxid = 'select accession, taxid from accession2taxid where accession in ("%s")'

    for acc_list in chunk_lists:
        filter = '","'.join(acc_list)
        data_annotation = cursor_refseq.execute(template_annotation % filter,).fetchall()
        data_taxid = cursor_taxid.execute(template_taxid % filter,).fetchall()
        accession2taxid = {}

        for row in data_taxid:
            accession2taxid[row[0]] = row[1]
        for annot in data_annotation:
            annot = list(annot)
            accession = annot[0]
            seq_length = annot[2]
            # AccessionVersion, TaxId, Title, Length
            description = re.sub("%s.[0-9]+ " % annot[0], "", annot[1])

            if accession in accession2taxid:
                taxid = accession2taxid[accession]
            else:
                taxid = accession2taxid_entrez(accession)
            hits.append([accession, taxid, description, seq_length])

    cursor.executemany(sql_template, hits)
    conn.commit()
    # index accession and taxid columns
    sql_index_1 = 'create index acc on refseq_hits (accession);'
    sql_index_2 = 'create index taxid on refseq_hits (taxid);'

    cursor.execute(sql_index_1)
    cursor.execute(sql_index_2)
    conn.commit()


def concatenate_core_orthogroups(fasta_files):
    out_name = 'msa.faa'

    # identification of all distinct fasta headers id (all unique taxons ids) in all fasta
    # storing records in all_seq_data (dico)
    taxons = []
    all_seq_data = {}
    for one_fasta in fasta_files:
        all_seq_data[one_fasta] = {}
        with open(one_fasta) as f:
            alignment = AlignIO.read(f, "fasta")
            for record in alignment:
                if record.id not in taxons:
                    taxons.append(record.id)
            all_seq_data[one_fasta][record.id] = record


    # building dictionnary of the form: dico[one_fasta][one_taxon] = sequence
    concat_data = {}

    start_stop_list = []
    start = 0
    stop = 0

    for one_fasta in fasta_files:
        start = stop + 1
        stop = start + len(all_seq_data[one_fasta][list(all_seq_data[one_fasta].keys())[0]]) - 1
        start_stop_list.append([start, stop])
        print(len(taxons))
        for taxon in taxons:
            # check if the considered taxon is present in the record
            if taxon not in all_seq_data[one_fasta]:
                 # if taxon absent, create SeqRecord object "-"*len(alignments): gap of the size of the alignment
                seq = Seq("-"*len(all_seq_data[one_fasta][list(all_seq_data[one_fasta].keys())[0]]))
                all_seq_data[one_fasta][taxon] = SeqRecord(seq, id=taxon)
            if taxon not in concat_data:
                concat_data[taxon] = all_seq_data[one_fasta][taxon]
            else:
                concat_data[taxon] += all_seq_data[one_fasta][taxon]

    # concatenating the alignments, writing to fasta file
    MSA = MultipleSeqAlignment([concat_data[i] for i in concat_data])
    with open(out_name, "w") as handle:
        AlignIO.write(MSA, handle, "fasta")

# Note : removed poor error handling
# (infinite recursion with infinite waiting time is 
# something goes wrong: better to crash fast than 
# a slow agonizing death)
def refseq_accession2fasta(accession_list):
    Entrez.email = "trestan.pillonel@unil.ch"
    handle = Entrez.efetch(db='protein', id=','.join(accession_list), rettype="fasta", retmode="text")
    records = [i for i in SeqIO.parse(handle, "fasta")]
    return records

def get_diamond_refseq_top_hits(databases_dir, phylum_filter, n_hits):
    conn = sqlite3.connect("orthology.db")
    cursor = conn.cursor()
    #sql1 = 'create table diamond_top_%s_nr_hits(orthogroup varchar, hit_accession, hit_sequence)' % ${params.refseq_diamond_BBH_phylogeny_top_n_hits}

    sql2 = 'attach "' + databases_dir + '/ncbi-taxonomy/linear_taxonomy.db" as linear_taxonomy'
    sql3 = 'attach "refseq_taxonomy.db" as refseq_taxonomy'
    sql4 = 'attach "diamond_refseq.db" as diamond_refseq'

    #cursor.execute(sql1,)
    cursor.execute(sql2,)
    cursor.execute(sql3,)
    cursor.execute(sql4,)
    conn.commit()

    #sql_template = 'insert into diamond_top_%s_nr_hits values (?,?,?)' % ${params.refseq_diamond_BBH_phylogeny_top_n_hits}

    sql_filtered_hits = """SELECT t1.orthogroup, t1.locus_tag, t3.sseqid, t3.hit_count 
       FROM locus_tag2orthogroup t1 INNER JOIN locus_tag2sequence_hash t2 ON t1.locus_tag=t2.locus_tag 
       INNER JOIN diamond_refseq.diamond_refseq t3 ON t2.sequence_hash=t3.qseqid 
       INNER JOIN refseq_taxonomy.refseq_hits t4 ON t3.sseqid=t4.accession
       INNER JOIN linear_taxonomy.ncbi_taxonomy t5 ON t4.taxid=t5.tax_id 
       WHERE t5.phylum not in ("%s")
       ORDER BY t1.orthogroup, t1.locus_tag, t3.hit_count;"""  % '","'.join(phylum_filter)
    filtered_hits = cursor.execute(sql_filtered_hits,).fetchall()

    # retrieve top hits excluding phylum of choice
    orthogroup2locus2top_hits = {}
    for row in filtered_hits:
        orthogroup = row[0]
        locus_tag = row[1]
        hit_id = row[2]

        if orthogroup not in orthogroup2locus2top_hits:
            orthogroup2locus2top_hits[orthogroup] = {}
        if locus_tag not in orthogroup2locus2top_hits[orthogroup]:
            orthogroup2locus2top_hits[orthogroup][locus_tag] = []
        if len(orthogroup2locus2top_hits[orthogroup][locus_tag]) < n_hits:
            orthogroup2locus2top_hits[orthogroup][locus_tag].append(hit_id)

    # retrieve aa sequences
    # why not select on orthogroups in orthogroup2locus2top_hits?
    # would be more efficient and spare some lines of code
    # as this is done in the loop below
    sql = """SELECT orthogroup, t1.locus_tag, sequence 
     FROM locus_tag2orthogroup t1 INNER JOIN locus_tag2sequence_hash t2 ON t1.locus_tag=t2.locus_tag 
     INNER JOIN sequence_hash2aa_sequence t3 ON t2.sequence_hash=t3.sequence_hash"""
    orthogroup_sequences = cursor.execute(sql,).fetchall()
    orthogroup2locus_and_sequence = {}
    for row in orthogroup_sequences:
        orthogroup = row[0]
        locus_tag = row[1]
        sequence = row[2]
        if orthogroup not in orthogroup2locus_and_sequence:
            orthogroup2locus_and_sequence[orthogroup] = []
        orthogroup2locus_and_sequence[orthogroup].append([locus_tag, sequence])

    # for each group, retrieve aa sequence from NCBI
    # write it to fasta file with orthogroup sequences
    for group in orthogroup2locus2top_hits:
        ortho_sequences = orthogroup2locus_and_sequence[group]
        refseq_sequence_records = []
        top_hits = [orthogroup2locus2top_hits[group][i] for i in orthogroup2locus2top_hits[group]]
        nr_top_hit = list(set([hit for hits_list in top_hits for hit in hits_list])) # flatten list... yeah, ugly.
        split_lists = chunks(nr_top_hit, 50)

        # multithreading here would be nice to run queries in parallel
        for one_list in split_lists:
            refseq_sequence_records += refseq_accession2fasta(one_list)

        if len(refseq_sequence_records) > 0:
            with open(group + "_nr_hits.faa", 'w') as f:
                for one_ortholog in ortho_sequences:
                    f.write(">%s\\n%s\\n" % (one_ortholog[0], one_ortholog[1]))
                for record in refseq_sequence_records:
                    name = record.name
                    f.write(">%s\\n%s\\n" % (name, str(record.seq)))

