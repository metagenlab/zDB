from Bio import Entrez, SeqIO
from Bio.SeqUtils import CheckSum
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from metagenlab_libs import db_utils

import setup_chlamdb


import pandas as pd
import itertools
import sys
import sqlite3
import urllib
import gzip
import re
import os
import datetime
import logging
import http.client

import mmap


# Heuristic to detect T3SS effectors inc proteins
# according to https://doi.org/10.1093/dnares/10.1.9
# - bi-lobed hydropathic domain, 55-65 aa long,
#   either 20-100 residues downstream the N-terminus or
#   100 residus from the C-terminus
# - no other hydrophic domain of more than 20 amino acid is present
# - the N-terminal portion doesn't contain any hydrophobic domain 
#    or signal peptide (to be determined)
# 
# Note : 
# - protein scale using the Kyte-Doolittle algorithm
#
# Input: fasta file
# Output: the set of proteins satisfying the conditions
# above

# Used to compute the hydrophobic plot
# CAVE : Should be an odd number
SLIDING_WINDOW = 7

N_TERMINUS_RANGE = (20, 100)
C_TERMINUS_LIMIT = 100 # must be within 100 aa from the C-terminal end
BILOBED_DOMAIN_MIN_SIZE = 55
BILOBED_DOMAIN_MAX_SIZE = 65


# if another hydrophobic domain bigger or equal to this
# is present, probably not an INC protein
HYDROPHOBIC_DOMAIN_EXCLUSION_SIZE = 20

# Amino-acid hydrophobicity values value, from
# Kyte J., Doolittle R.F. 
# J. Mol. Biol. 157:105-132(1982).
AA_HYDROPHOBICITY_SCALE_VALUES = {
    'A' : 1.8,
    'R' : -4.5,
    'N' : -3.5,
    'D' : -3.5,
    'C' : 2.5,
    'Q' : -3.5,
    'E' : -3.5,
    'G' : -0.4,
    'H' : -3.2,
    'I' : 4.5,
    'L' : 3.8,
    'K' : -3.9,
    'M' : 1.9,
    'F' : 2.8,
    'P' : -1.6,
    'S' : -0.8,
    'T' : -0.7,
    'W' : -0.9,
    'Y' : -1.3,
    'V' : 4.2
}


# Note: the criteria are a bit loosened as candidates will
# necessitate a manual control of the hydropathy plot anyway
#
# TODO : add detection of the different types of effectors according to the paper
# Type I, II and III and modify the manual handpicking accordingly
def T3SS_inc_proteins_detection(fasta_file, out_file):
    T3SS_hydropathy_values = []

    # because of the sliding window algorithm, we have results for the following
    # positions [index_shift;index_shift + 1... ; len-index_shift-1; len-index_shift]
    index_shift = (SLIDING_WINDOW+1)/2

    for record in SeqIO.parse(fasta_file, "fasta"):
        analysis = ProteinAnalysis(str(record.seq))
        values = analysis.protein_scale(AA_HYDROPHOBICITY_SCALE_VALUES, SLIDING_WINDOW)

        hydropathy_domains = []
        i = 0
        while i<len(values):
            # not hydrophic, skip
            if values[i] < 0:
                i += 1
                continue
            
            j = i+1
            # Note : may be too stringent, what if a single value
            # in the middle of an otherwise perfect sequence is
            # too low. Should we keep it or tolerate short hydrophilic window
            # in the middle of an overall hydrophobic bilobed domain ?
            while j<len(values) and values[j]>=0:
                j+=1

            # not interesting to take into account if size
            # lower than this --> wouldn't change anything
            if j-i>=HYDROPHOBIC_DOMAIN_EXCLUSION_SIZE:
                hydropathy_domains.append((i, j-1))
            i = j+1

        C_terminus_start_limit = len(values)-C_TERMINUS_LIMIT
        n_domains = 0
        has_other_hydrophobic_domain = False
        bilobed_domain = []
        for start, end in hydropathy_domains:
            length = end-start+1
            if length<BILOBED_DOMAIN_MIN_SIZE or length>BILOBED_DOMAIN_MAX_SIZE:
                has_other_hydrophobic_domain = True
                break

            # see if end or beginning of the domain is within bounds
            # Note: this is less stringent than the initial conditions
            if ((start+index_shift>=N_TERMINUS_RANGE[0] and start+index_shift<=N_TERMINUS_RANGE[1]) \
                    or (end+index_shift>=N_TERMINUS_RANGE[0] and end+index_shift<=N_TERMINUS_RANGE[1]) \
                    or end>=C_terminus_start_limit+index_shift):
                n_domains += 1
                bilobed_domain = (start, end)
            else:
                has_other_hydrophobic_domain = True
                break

        if n_domains==1 and not has_other_hydrophobic_domain:
            T3SS_hydropathy_values.append([values, bilobed_domain, record.id])
    
    hydrophobic_plot_file = open(out_file, "w")
    for values, domain, record_id in T3SS_hydropathy_values:
        string_val = [str(i) for i in values]
        hydrophobic_plot_file.write(f">{record_id}\n")
        hydrophobic_plot_file.write(f">{domain[0]}-{domain[1]}\n")
        hydrophobic_plot_file.write("\n".join(string_val) + "\n")

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]

def get_idmapping_crossreferences(databases_dir, table):
    conn = sqlite3.connect(databases_dir + "/uniprot/idmapping/idmapping.db")
    cursor = conn.cursor()
    o = open("idmapping_crossreferences.tab", "w")
    sql = """SELECT uniprokb_accession, db_name, accession 
        FROM uniparc2uniprotkb t1 INNER JOIN uniprotkb_cross_references t2 ON t1.uniprotkb_id=t2.uniprotkb_id 
        INNER JOIN database t3 ON t2.db_id=t3.db_id 
        WHERE uniparc_accession=?;"""

    with open(table, 'r') as f:
        f.readline()
        for row in f:
            data = row.rstrip().split("\t")
            uniparc_accession = data[2]
            checksum = data[0]
            cursor.execute(sql, [uniparc_accession])
            crossref_list = cursor.fetchall()
            for crossref in crossref_list:
                uniprot_accession = crossref[0]
                db_name = crossref[1]
                db_accession = crossref[2]
                o.write("%s\t%s\t%s\t%s\n" % ( checksum,
                                               uniprot_accession,
                                               db_name,
                                               db_accession))

def get_uniparc_crossreferences(databases_dir, table):
    conn = sqlite3.connect(databases_dir + "/uniprot/uniparc/uniparc.db")
    cursor = conn.cursor()
    o = open("uniparc_crossreferences.tab", "w")

    sql = """SELECT db_name, accession, status 
        FROM uniparc_cross_references t1 INNER JOIN crossref_databases t2 ON t1.db_id=t2.db_id 
        WHERE uniparc_id=? AND db_name NOT IN ("SEED", "PATRIC", "EPO", "JPO", "KIPO", "USPTO");"""

    with open(table, 'r') as f:
        # skip header
        f.readline()
        for row in f:
            data = row.rstrip().split("\t")
            uniparc_id = str(data[1])
            checksum = data[0]
            cursor.execute(sql, [uniparc_id])
            crossref_list = cursor.fetchall()
            for crossref in crossref_list:
                db_name = crossref[0]
                db_accession = crossref[1]
                entry_status = crossref[2]
                o.write("%s\t%s\t%s\t%s\n" % ( checksum,
                                          db_name,
                                          db_accession,
                                          entry_status))


def get_oma_mapping(databases_dir, fasta_file):
    conn = sqlite3.connect(databases_dir + "/oma/oma.db")
    cursor = conn.cursor()
    oma_map = open('oma_mapping.tab', 'w')
    no_oma_mapping = open('no_oma_mapping.faa', 'w')
    oma_map.write("locus_tag\toma_id\n")
    records = SeqIO.parse(fasta_file, "fasta")
    no_oma_mapping_records = []

    # TODO : parallelize database requests or make a single big request
    #  using IN with all the accessions
    for record in records:
        sql = 'select accession from hash_table where sequence_hash=?'
        cursor.execute(sql, (CheckSum.seguid(record.seq),))
        hits = cursor.fetchall()
        if len(hits) == 0:
            no_oma_mapping_records.append(record)
        for hit in hits:
            oma_map.write("%s\t%s\n" % (record.id, hit[0]))
    SeqIO.write(no_oma_mapping_records, no_oma_mapping, "fasta")

def pmid2abstract_info(pmid_list):
    from Bio import Medline

    # make sure that pmid are strings
    pmid_list = [str(i) for i in pmid_list]

    try:
        handle = Entrez.efetch(db="pubmed", id=','.join(pmid_list), rettype="medline", retmode="text")
        records = Medline.parse(handle)
    except:
        print("FAIL:", pmid)
        return None

    pmid2data = {}
    for record in records:
        try:
            pmid = record["PMID"]
        except:
            print(record)
            #{'id:': ['696885 Error occurred: PMID 28696885 is a duplicate of PMID 17633143']}
            if 'duplicate' in record['id:']:
                duplicate = record['id:'].split(' ')[0]
                correct = record['id:'].split(' ')[-1]
                print("removing duplicated PMID... %s --> %s" % (duplicate, correct))
                # remove duplicate from list
                pmid_list.remove(duplicate)
                return pmid2abstract_info(pmid_list)

        pmid2data[pmid] = {}
        pmid2data[pmid]["title"] = record.get("TI", "?")
        pmid2data[pmid]["authors"] = record.get("AU", "?")
        pmid2data[pmid]["source"] = record.get("SO", "?")
        pmid2data[pmid]["abstract"] = record.get("AB", "?")
        pmid2data[pmid]["pmid"] = pmid
    return pmid2data

def get_PMID_data():
    Entrez.email = "trestan.pillonel@chuv.ch"
    conn = sqlite3.connect("string_mapping_PMID.db")
    cursor = conn.cursor()
    sql = 'create table if not exists hash2pmid (hash binary, pmid INTEGER)'
    sql2 = 'create table if not exists pmid2data (pmid INTEGER, title TEXT, authors TEXT, source TEXT, abstract TEXT)'
    cursor.execute(sql,)
    cursor.execute(sql2,)

    # get PMID nr list and load PMID data into sqldb
    pmid_nr_list = []
    sql_template = 'insert into hash2pmid values (?, ?)'
    with open("string_mapping_PMID.tab", "r") as f:
        n = 0
        for row in f:
            data = row.rstrip().split("\t")
            if data[1] != 'None':
                n+=1
                cursor.execute(sql_template, data)
                if data[1] not in pmid_nr_list:
                    pmid_nr_list.append(data[1])
            if n % 1000 == 0:
                print(n, 'hash2pmid ---- insert ----')
                conn.commit()

    pmid_chunks = [i for i in chunks(pmid_nr_list, 50)]

    # get PMID data and load into sqldb
    sql_template = 'insert into pmid2data values (?, ?, ?, ?, ?)'
    for n, chunk in enumerate(pmid_chunks):
        print("pmid2data -- chunk %s / %s" % (n, len(pmid_chunks)))
        pmid2data = pmid2abstract_info(chunk)
        for pmid in pmid2data:
            cursor.execute(sql_template, (pmid, pmid2data[pmid]["title"], str(pmid2data[pmid]["authors"]), pmid2data[pmid]["source"], pmid2data[pmid]["abstract"]))
        if n % 10 == 0:
            conn.commit()
    conn.commit()

# Removed error handling possibly leading
# to endless loop. A quick death is usually better.
# TODO: we should improve the error handling rather than 
# remove it (with eg a max number of retry) 

def uniprot_accession2score(uniprot_accession_list):
    # https://www.uniprot.org/uniprot/?query=id:V8TQN7+OR+id:V8TR74&format=xml
    link = "http://www.uniprot.org/uniprot/?query=id:%s&columns=id,annotation%%20score&format=tab" % ("+OR+id:".join(uniprot_accession_list))

    page = urllib.request.urlopen(link)
    data = page.read().decode('utf-8').split('\n')
    rows = [i.rstrip().split('\t') for i in data]
    unirpot2score = {}
    for row in rows:
        if len(row) > 0:
            if row[0] == 'Entry':
                continue
            elif len(row)<2:
                continue
            else:
                unirpot2score[row[0]] = row[1]
    return unirpot2score


def get_uniprot_data(databases_dir, table):
    conn = sqlite3.connect(databases_dir + "/uniprot/idmapping/uniprot_sprot_trembl.db")
    cursor = conn.cursor()
    uniprot_table = open(table, 'r')
    uniprot_data = open('uniprot_data.tab', 'w')
    uniprot_data.write("uniprot_accession\tuniprot_score\tuniprot_status\tproteome\tcomment_function\tec_number\tcomment_subunit\tgene\trecommendedName_fullName\tproteinExistence\tdevelopmentalstage\tcomment_similarity\tcomment_catalyticactivity\tcomment_pathway\tkeywords\n")

    uniprot_accession_list = [row.rstrip().split("\t")[1].split(".")[0] for row in uniprot_table if row.rstrip().split("\t")[1] != 'uniprot_accession']
    uniprot_accession_chunks = chunks(uniprot_accession_list, 300)

    sql_uniprot_annot = 'SELECT * FROM uniprot_annotation WHERE uniprot_accession IN ("%s");'

    # parallelizing those requests could be interesting if this function
    # took too much time
    for one_chunk in uniprot_accession_chunks:
        uniprot2score = uniprot_accession2score(one_chunk)
        filter = '","'.join(one_chunk)
        uniprot_annot_data = cursor.execute(sql_uniprot_annot % filter,).fetchall()

        for uniprot_annotation in uniprot_annot_data:
            uniprot_accession = uniprot_annotation[0]
            if uniprot_accession in uniprot2score:
                uniprot_score = uniprot2score[uniprot_accession]
            else:
                uniprot_score = 0
            comment_function = uniprot_annotation[1]
            ec_number = uniprot_annotation[2]
            comment_similarity = uniprot_annotation[3]
            comment_catalyticactivity = uniprot_annotation[4]
            comment_pathway = uniprot_annotation[5]
            keywords = uniprot_annotation[6]
            comment_subunit = uniprot_annotation[7]
            gene = uniprot_annotation[8]
            recommendedName_fullName = uniprot_annotation[9]
            proteinExistence = uniprot_annotation[10]
            developmentalstage = uniprot_annotation[11]
            proteome = uniprot_annotation[12]
            reviewed = uniprot_annotation[14]

            uniprot_data.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                    % (uniprot_accession, uniprot_score, reviewed, proteome, comment_function, ec_number,
                       comment_subunit, gene, recommendedName_fullName, proteinExistence,
                       developmentalstage, comment_similarity, comment_catalyticactivity,
                       comment_pathway, keywords))

def get_uniparc_mapping(databases_dir, fasta_file):
    conn = sqlite3.connect(databases_dir + "/uniprot/uniparc/uniparc.db")
    cursor = conn.cursor()

    uniparc_map = open('uniparc_mapping.tab', 'w')
    uniprot_map = open('uniprot_mapping.tab', 'w')
    no_uniprot_mapping = open('no_uniprot_mapping.faa', 'w')
    no_uniparc_mapping = open('no_uniparc_mapping.faa', 'w')
    uniparc_mapping_faa = open('uniparc_mapping.faa', 'w')

    uniparc_map.write("locus_tag\tuniparc_id\tuniparc_accession\tstatus\n")
    uniprot_map.write("locus_tag\tuniprot_accession\ttaxon_id\tdescription\n")

    records = SeqIO.parse(fasta_file, "fasta")
    no_mapping_uniprot_records = []
    no_mapping_uniparc_records = []
    mapping_uniparc_records = []

    for record in records:
        match = False
        sql = """SELECT t1.uniparc_id, uniparc_accession, accession,taxon_id, description, db_name, status
            FROM uniparc_accession t1 INNER JOIN uniparc_cross_references t2 ON t1.uniparc_id=t2.uniparc_id
            INNER JOIN crossref_databases t3 ON t2.db_id=t3.db_id 
            WHERE sequence_hash=?"""
        cursor.execute(sql, (CheckSum.seguid(record.seq),))
        hits = cursor.fetchall()
        if len(hits) == 0:
            no_mapping_uniparc_records.append(record)
            no_mapping_uniprot_records.append(record)
        else:
            # All rows have the same UP 
            # accession, but we have to check that there is at least one active cross-reference.
            # If all entries are dead, it means that the sequence was removed from 
            # all cross-referenced databases. In that case, we consider it as a "new"
            # sequence and annotate it "de novo" with interproscan.
            mapping_uniparc_records.append(record)
            all_status = [i[6] for i in hits]
            if 1 in all_status:
                status = 'active'
            else:
                status = 'dead'
            uniparc_map.write(f"{record.id}\t{hits[0][0]}\t{hits[0][1]}\t{status}\n")
            for uniprot_hit in hits:
                if uniprot_hit[5] in ["UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL"] and uniprot_hit[6] == 1:
                    match = True
                    uniprot_map.write("%s\t%s\t%s\t%s\t%s\n" % (record.id,
                                                                 uniprot_hit[2],
                                                                 uniprot_hit[3],
                                                                 uniprot_hit[4],
                                                                 uniprot_hit[5]))
            if not match:
                no_mapping_uniprot_records.append(record)
    SeqIO.write(no_mapping_uniprot_records, no_uniprot_mapping, "fasta")
    SeqIO.write(no_mapping_uniparc_records, no_uniparc_mapping, "fasta")
    SeqIO.write(mapping_uniparc_records, uniparc_mapping_faa, "fasta")

def refseq_locus_mapping(accessions_list):
    o = open("refseq_corresp.tab", "w")
    for accession in accessions_list:
        with gzip.open(accession, "rt") as handle:
            records = SeqIO.parse(handle, "genbank")
            for record in records:
                for feature in record.features:
                    if 'protein_id' in feature.qualifiers and 'old_locus_tag' in feature.qualifiers:
                        refseq_locus_tag = feature.qualifiers["locus_tag"][0]
                        protein_id = feature.qualifiers["protein_id"][0]
                        old_locus_tag = feature.qualifiers["old_locus_tag"][0]
                        o.write(f"{old_locus_tag}\t{refseq_locus_tag}\t{protein_id}\n")


def merge_plasmids(plasmids_records):
    """
    Same idea as for merge_gbk: some plasmids may be composed of several contigs.
    The identification is based on the name in the /plasmid qualifier in the feature:
    if the name is the same, then the contigs are grouped with a 200*N assembly gap and the 
    informations from the first contig is kept.
    """

    hsh_name_to_records = {}
    untouched_records = []

    # group the records in fonction of the plasmid name
    for record in plasmids_records:
        found_name = False
        for feature in record.features:
            if feature.type != "source":
                continue
            if "plasmid" not in feature.qualifiers:
                continue
            names = feature.qualifiers["plasmid"]
            if len(names) == 0:
                continue
            name = names[0]
            lst = hsh_name_to_records.get(name, [])
            lst.append(record)
            hsh_name_to_records[name] = lst
            found_name = True
            break

        if not found_name:
            untouched_records.append(record)

    modified_records = []
    for name, records in hsh_name_to_records.items():
        modified_records.append(merge_gbk(records, plasmids=True))

    return untouched_records + modified_records


def merge_gbk(gbk_records, filter_size=0, gi=False, plasmids=False):
    '''
    merge multiple contigs into a single DNA molecule with 200*N between contigs
    keep source description from the first record
    remove contigs smaller than <filter_size>

    For plasmids: keep the first accession, as having several entries with the 
    same accession would profoundly hurt BioSQL.

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

    accession_index = 0 if plasmids else -1
    try:
        merged_rec.id = gbk_records[0].annotations["accessions"][accession_index]
    except KeyError:
        merged_rec.id = gbk_records[0].id

    if gi:
        merged_rec.annotations["gi"] = gi

    merged_rec.description = "%s" % gbk_records[0].annotations["organism"]
    merged_rec.annotations = gbk_records[0].annotations
    try:
        merged_rec.name = gbk_records[0].annotations["accessions"][accession_index]
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
        if setup_chlamdb.is_plasmid(record):
            plasmid_record_list.append(record)
        else:
            chromosome_record_list.append(record)
    return chromosome_record_list, plasmid_record_list


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


def orthogroups_to_fasta(genomes_list):
    fasta_list = genomes_list.split(' ')

    sequence_data = {}
    for fasta_file in fasta_list:
        sequence_data.update(SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta")))

      # write fasta
    with open("Orthogroups.txt") as f:
        all_grp = [i for i in f]
        for n, line in enumerate(all_grp):
            groups = line.rstrip().split(' ')
            group_name = groups[0][:-1]
            groups = groups[1:]
            if len(groups)>1:
                new_fasta = [sequence_data[i] for i in groups]
                out_path = "%s.faa" % group_name
                out_handle = open(out_path, "w")
                SeqIO.write(new_fasta, out_handle, "fasta")


def check_gbk(gbff_files):
    reannotation_list = []

    # count the number of identical source names
    source2count = {}
    accession2count = {}
    for genome_number, gbff_file in enumerate(gbff_files):

        records = list(SeqIO.parse(gbff_file, "genbank"))
        for record in records:
            n_missing, total = count_missing_locus_tags(record)
            if n_missing > 0:
                print ('Warning: %s/%s missing locus tag for record %s' % (n_missing, total, record.name))
        chromosome, plasmids = filter_plasmid(records)
        cleaned_records = []
        plasmid_reannot = False
        chromosome_reannot = False

        for n_plasmid, plasmid in enumerate(plasmids):
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
            '''

            annot_check = is_annotated(chromosome[0])
            if annot_check:
                if len(chromosome) > 1:
                    merged_record = merge_gbk(chromosome)
                else:
                    merged_record = chromosome[0]
                merged_record.description = clean_description(merged_record.description)
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


def check_organism_uniqueness(gbk_lst):
    """
    As BioSQL uses the organism entry of the records to assign
    a taxid when it is not allowed to access the ncbi online,
    we need to ensure that the same organism name is not used
    twice in the genomes to prevent bioentries from different genbank
    files to be entered under the same taxon_id (as chlamdb uses the taxon_id to 
    differentiate between genomes).

    Display an error message and stop the pipeline if this arises
    """
    organisms = dict()

    for gbk_file in gbk_lst:
        curr_organism = None
        for record in SeqIO.parse(gbk_file, "genbank"):
            if not "organism" in record.annotations:
                raise Exception(f"No organism for file {gbk_file}")
            
            organism = record.annotations["organism"]
            if curr_organism is None:
                curr_organism = organism
            elif curr_organism != organism:
                raise Exception(f"Two different organism in {gbk_file}: {curr_organism}/{organism}")

        gbk_lst = organisms.setdefault(organism, [])
        gbk_lst.append(gbk_file)

    for organism, file_lst in organisms.items():
        if len(file_lst) > 1:
            genbanks = ",".join(file_lst)
            raise Exception(f"Genbank files {genbanks} have the same organism: {organism}")


def check_names(record, common_names):
    common_name = record.annotations.get("source", None)
    if (common_name is None) or (common_name in common_names):
        # such a common name has already been 
        # encountered: change it to the unique scientific
        # name
        # NOTE: We know that record.annotations contains "organism"
        record.annotations["source"] = record.annotations["organism"]
    return record


def filter_out_unannotated(gbk_files):
    """
    NOTE: assumes that the check_organism_uniqueness has been called before (i.e.
    unique scientific organism name and only one organism per file).

    This function does two things:
    ** filter out unannotated contigs (TODO: check if this is really necessary)
    ** for genbank files that have the same organism common name, change it to their
       scientific name (BioSQL would ortherwise assign them the same taxon_id,
       I know, this is stupid)
    """
    common_names = set()
    for gbk_file in gbk_files:
        records = SeqIO.parse(gbk_file, "genbank")
        result_file = gbk_file.replace(".gbk", "_filtered.gbk")
        to_keep = (check_names(rec, common_names) for rec in records if is_annotated(rec))
        first_rec = next(to_keep)
        common_names.add(first_rec.annotations["source"])
        SeqIO.write(itertools.chain([first_rec], to_keep), result_file, "genbank")


list_locus_tag_nr=[]
list_locus_tag=[]

def lt_uniqueness_gbk(gbk_list):
 for gbk_file in gbk_list:
  records = SeqIO.parse(gbk_file, 'genbank')  #from BioPython, object of class SeqRecord
  genome = os.path.basename(gbk_file).split('.')[0]

  for record in records:
   record_id=record.id
   for feature in record.features:
     if (feature.type == 'CDS' and 'pseudo' not in feature.qualifiers and 'pseudogene' not in feature.qualifiers):
       if "locus_tag" in feature.qualifiers:
          locus_tag = feature.qualifiers["locus_tag"][0]
          list_locus_tag.append(locus_tag)
          count_locus_tag= len(list_locus_tag)
          if locus_tag not in list_locus_tag_nr:
                list_locus_tag_nr.append(locus_tag)
                count_nr_locus_tag= len(list_locus_tag_nr)
          if locus_tag=="":
             raise Exception(f"Empty locus_tag in {record_id} of {genome}")
       else:
             raise Exception(f"No locus_tag in {record_id} of {genome}")

 if not count_locus_tag==count_nr_locus_tag:
  raise Exception(f"Non-unique locus_tag")

def checked_unique_lt (gbk_files):
 for gbk_file in gbk_files:
  records = SeqIO.parse(gbk_file, 'genbank')
  new_gbk = gbk_file.replace(".gbk", "_lt_checked.gbk") 
  string= "list"
  list= string.join({new_gbk})
  list=[]
  for record in records:
   list.append(record)
  pass
  SeqIO.write( list, new_gbk, 'genbank')





LOCUS_id_list=[]
LOCUS_id_list_nr=[]

def LOCUS_uniqueness_gbk(gbk_list):
 for gbk_file in gbk_list:
  records = SeqIO.parse(gbk_file, 'genbank')  #from BioPython, object of class SeqRecord
  count_feature_in_record=0

  for record in records:
   LOCUS_id=record.name
   for feature in record.features:
     if (feature.type == 'CDS' and 'pseudo' not in feature.qualifiers and 'pseudogene' not in feature.qualifiers):
       count_feature_in_record=count_feature_in_record + 1
   if count_feature_in_record>=1:
    LOCUS_id_list.append(LOCUS_id)
    count_LOCUS_id=len(LOCUS_id_list)
    if LOCUS_id not in LOCUS_id_list_nr:
       LOCUS_id_list_nr.append(LOCUS_id)
       count_nr_LOCUS_id= len(LOCUS_id_list_nr)
 if not  count_nr_LOCUS_id==count_LOCUS_id:
   raise Exception(f"Non-unique LOCUS identifier")


def checked_unique_L (gbk_files):
 for gbk_file in gbk_files:
  records = SeqIO.parse(gbk_file, 'genbank')
  new_gbk = gbk_file.replace(".gbk", "_lt_L_checked.gbk") 
  string= "list"
  list= string.join({new_gbk})
  list=[]
  for record in records:
   list.append(record)
  pass
  SeqIO.write( list, new_gbk, 'genbank')



def string_id2pubmed_id_list(accession):
    link = 'http://string-db.org/api/tsv/abstractsList?identifiers=%s' % accession
    try:
        data = urllib.request.urlopen(link).read().rstrip().decode('utf-8').split('\n')[1:]
    except urllib.URLError:
        return False
    pid_list = [row.split(':')[1] for row in data]
    return pid_list


def get_string_PMID_mapping(string_map):
    o = open("string_mapping_PMID.tab", "w")
    with open(string_map, 'r') as f:
        # skip header
        f.readline()
        for row in f:
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
        for hit in hits:
            pdb_map.write("%s\t%s\n" % (record.id, hit[0]))
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
        for hit in hits:
            tcdb_map.write("%s\t%s\n" % (record.id, hit[0]))

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


def convert_gbk_to_faa_SEQ(gbf_file, edited_gbf):
    records = SeqIO.parse(gbf_file, 'genbank')  #from BioPython, object of class SeqRecord
    edited_records = open(edited_gbf, 'w')
    genome = (gbf_file).split('.')[0]


    for record in records:
        description=record.description
        record_name=record.name
        all_trans_seq=""

        for feature in record.features:
           if "translation" in feature.qualifiers:
             translation_seq = feature.qualifiers["translation"][0]
             translation_seq_str=''.join(translation_seq)
             all_trans_seq= all_trans_seq + translation_seq_str
        edited_records.write(">%s %s\n%s\n" % (record_name, description, all_trans_seq))



def convert_gbk_to_fna_SEQ (gbf_file, fna_contigs):
    records = SeqIO.parse(gbf_file, 'genbank')  #from BioPython, object of class SeqRecord
    edited_records = open(fna_contigs, 'w')

    for record in records:
       edited_records.write(">%s %s\n%s\n" % (record.name, record.description, record.seq))




def convert_gbk_to_ffn_seq (gbf_file, edited_gbf):
    records = SeqIO.parse(gbf_file, 'genbank')  
    edited_records = open(edited_gbf, 'w')

    for record in records:
        protein2count = {}
        record_id = record.id
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
                                                     record.description, feature.location.extract(record).seq))
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

# all faa files are merged into fasta_file
def get_nr_sequences(fasta_file, genomes_list):
    locus2genome = {}
    for fasta in genomes_list:
        genome = os.path.basename(fasta).split('.')[0]
        for seq in SeqIO.parse(fasta, "fasta"):
            locus2genome[seq.name] = genome
    nr_fasta = open('nr.faa', 'w')
    nr_mapping = open('nr_mapping.tab', 'w')

    hsh_checksum_list = {}

    records = SeqIO.parse(fasta_file, "fasta")
    updated_records = []

    for record in records:

        # NOTE: the case is important for crc64, need to check whether it
        # is necessary to make all entries lower/upper case to ensure consistency.
        checksum = CheckSum.crc64(record.seq)
        nr_mapping.write("%s\t%s\t%s\n" % (record.id,
                                          checksum,
                                          locus2genome[record.id]))
        if checksum not in hsh_checksum_list:
            hsh_checksum_list[checksum] = [record]
            record.id = checksum
            record.name = ""
            updated_records.append(record)
        else:
            # NOTE: having same hash does not mean that the sequences are identical: as
            # the hash space is smaller than the sequence space, it means that collision
            # are unavoidable (but not probable) and record with same hashes should be compared
            # https://www.uniprot.org/help/uniparc (sequence comparison)
            #
            # the list of records having the same checksum, but potentially, 
            # different sequences -> compare them: python does so 
            # comparing the sequences as strings, assuming a similar alphabet
            lst_records = hsh_checksum_list[checksum]
            sequence = record.seq
            has_identical = False
            for prev_record in lst_records:
                if prev_record.seq == sequence:
                    has_identical = True
                    break
            if not has_identical:
                lst_records.append(record)
                record.id = checksum + "-" + len(lst_records)
                record.name = ""
                updated_records.append(record)

    SeqIO.write(updated_records, nr_fasta, "fasta")



# This function parses the result of orthofinder/orthoMCL, 
# retrieves the groups of orthologs detected by the tools
# and returns the core single copy orthologs (i.e. the 
# ortholog present in all samples)
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

    df =df.apply(pd.to_numeric, args=('coerce',))

    n_genomes = len(set(locus2genome.values()))
    n_minimum_genomes = n_genomes-n_missing
    freq_missing = (n_genomes-float(n_missing))/n_genomes
    limit = freq_missing*n_genomes

    groups_with_paralogs = df[(df > 1).sum(axis=1) > 0].index
    df = df.drop(groups_with_paralogs)

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


def concatenate_core_orthogroups(fasta_files):
    out_name = 'msa.faa'

    # identification of all distinct fasta headers id (all unique taxons ids) in all fasta
    # storing records in all_seq_data (dico)
    taxons = []
    all_seq_data = {}
    for one_fasta in fasta_files:
        all_seq_data[one_fasta] = {}
        for record in AlignIO.read(one_fasta, "fasta"):
            if record.id not in taxons:
                taxons.append(record.id)
            all_seq_data[one_fasta][record.id] = record

    # building dictionnary of the form: dico[one_fasta][one_taxon] = sequence
    concat_data = {}
    for one_fasta in fasta_files:
        for taxon in taxons:
            # since core orthogroups, all taxids should be present everywhere
            assert(taxon in all_seq_data[one_fasta])
            if taxon not in concat_data:
                concat_data[taxon] = all_seq_data[one_fasta][taxon]
            else:
                concat_data[taxon] += all_seq_data[one_fasta][taxon]

    # concatenating the alignments, writing to fasta file
    MSA = MultipleSeqAlignment([concat_data[i] for i in concat_data])
    with open(out_name, "w") as handle:
        AlignIO.write(MSA, handle, "fasta")

# TODO: we should improve the error handling rather than remove it 
# because errors can occur (kind of randomly) with Entez (those errors
# are not reproducible, this is why I did not care dealing with the 
# infinite call)
# NOTE: redundant code, kept for now until chlamdb is fully functional with 
# the new version
def refseq_accession2fasta(accession_list):
    Entrez.email = "trestan.pillonel@chuv.ch"
    handle = Entrez.efetch(db='protein', id=','.join(accession_list), rettype="fasta", retmode="text")
    records = [i for i in SeqIO.parse(handle, "fasta")]
    return records


def get_uniprot_goa_mapping(database_dir, uniprot_acc_list):
    conn = sqlite3.connect(database_dir + "/goa/goa_uniprot.db")
    cursor = conn.cursor()
    o = open("goa_uniprot_exact_annotations.tab", "w")

    sql = """SELECT GO_id, reference, evidence_code, category 
        FROM goa_table 
        WHERE uniprotkb_accession IN (\"{}\");"""

    accessions_base = [i.split(".")[0].strip() for i in uniprot_acc_list]
    query = sql.format("\",\"".join(accessions_base))
    cursor.execute(sql.format(",".join(accessions_base)),)
    go_list = cursor.fetchall()
    for go in go_list:
        GO_id = go[0]
        reference = go[1]
        evidence_code = go[2]
        category = go[3]
        o.write(f"{accession_base}\t{GO_id}\t{reference}\t{evidence_code}\t{category}\n")

