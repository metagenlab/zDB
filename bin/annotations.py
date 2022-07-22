from Bio import Entrez, SeqIO
from Bio.SeqUtils import CheckSum
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import logging
import sys
import sqlite3
import urllib
import gzip
import re
import os
import datetime
import logging
import http.client
import pandas as pd
from ftplib import FTP

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
            # NOTE : may be too stringent, what if a single value
            # in the middle of an otherwise perfect sequence is
            # too low? Should we keep it or tolerate short hydrophilic window
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


# def T3SS_inc_protein_interproscan_based(interpro_output_files, out_file):

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
    # https://rest.uniprot.org/uniprotkb/?query=id:V8TQN7+OR+id:V8TR74&format=xml
    link = "https://rest.uniprot.org/uniprotkb/search?query=accession:%s&fields=accession,annotation_score&format=tsv" % ("+OR+accession:".join(uniprot_accession_list))

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
    # uniprot_id  uniprot_accession  comment_function                                              
    # ec_number  comment_similarity                        
    # comment_catalyticactivity  comment_pathway  keywords  comment_subunit  gene  recommendedName_fullName  proteinExistence        developmentalstage  proteome     annotation_score  reviewed
    for one_chunk in uniprot_accession_chunks:
        uniprot2score = uniprot_accession2score(one_chunk)
        filter = '","'.join(one_chunk)
        uniprot_annot_data = pd.read_sql(sql_uniprot_annot % filter, conn)

        for n,uniprot_annotation in uniprot_annot_data.iterrows():
            uniprot_accession = uniprot_annotation.uniprot_accession
            if uniprot_accession in uniprot2score:
                uniprot_score = uniprot2score[uniprot_accession]
            else:
                uniprot_score = 0
            comment_function = uniprot_annotation.comment_function
            ec_number = uniprot_annotation.ec_number
            comment_similarity = uniprot_annotation.comment_similarity
            comment_catalyticactivity = uniprot_annotation.comment_catalyticactivity
            comment_pathway = uniprot_annotation.comment_pathway
            keywords = uniprot_annotation.keywords
            comment_subunit = uniprot_annotation.comment_subunit
            gene = uniprot_annotation.gene
            recommendedName_fullName = uniprot_annotation.recommendedName_fullName
            proteinExistence = uniprot_annotation.proteinExistence
            developmentalstage = uniprot_annotation.developmentalstage
            proteome = uniprot_annotation.proteome
            reviewed = uniprot_annotation.reviewed

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
            # All rows have the same UP accession, but we have to check that there 
            # is at least one active cross-reference.
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

def download_assembly_refseq(accession):
    Entrez.email = "trestan.pillonel@chuv.ch"
    Entrez.api_key = "719f6e482d4cdfa315f8d525843c02659408"
    if accession=="":
        return

    handle1 = Entrez.esearch(db="assembly", term=accession)
    record1 = Entrez.read(handle1)
    ncbi_id = record1['IdList'][-1]
    handle_assembly = Entrez.esummary(db="assembly", id=ncbi_id)
    assembly_record = Entrez.read(handle_assembly, validate=False)

    # only consider annotated genbank => get corresponding refseq assembly
    if 'genbank_has_annotation' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]["PropertyList"]:
        ftp_path = re.findall('<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
        ftp=FTP('ftp.ncbi.nih.gov')
        ftp.login("anonymous","trestan.pillonel@unil.ch")
        ftp.cwd(ftp_path)
        filelist=ftp.nlst()
        filelist = [i for i in filelist if 'genomic.gbff.gz' in i]
        for file in filelist:
            ftp.retrbinary("RETR "+file, open(file, "wb").write)
    else:
        print("no genbank annotation for ${accession}")

def download_assembly(accession):
    Entrez.email = "trestan.pillonel@chuv.ch"
    Entrez.api_key = "719f6e482d4cdfa315f8d525843c02659408"

    handle1 = Entrez.esearch(db="assembly", term=accession)
    record1 = Entrez.read(handle1)

    ncbi_id = record1['IdList'][-1]
    handle_assembly = Entrez.esummary(db="assembly", id=ncbi_id)
    assembly_record = Entrez.read(handle_assembly, validate=False)

    if 'genbank_has_annotation' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]["PropertyList"]:
        ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
    elif 'refseq_has_annotation' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]["PropertyList"]:
      ftp_path = re.findall('<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
    else:
        raise("${accession} assembly not annotated! --- exit ---")
    ftp=FTP('ftp.ncbi.nih.gov')
    ftp.login("anonymous","trestan.pillonel@unil.ch")
    ftp.cwd(ftp_path)
    filelist=ftp.nlst()
    filelist = [i for i in filelist if 'genomic.gbff.gz' in i]
    for file in filelist:
        ftp.retrbinary("RETR " + file, open(file, "wb").write)


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

# See documentation of GenBank record format for more informations
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
    to_keep = [rec for rec in SeqIO.parse(gbk_file, "genbank") if is_annotated(rec)]
    SeqIO.write(to_keep, result_file, "genbank")


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

# filter out small sequences and ambiguous amino-acids
def filter_sequences(fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    processed_records = []
    for record in records:
        if len(record.seq) >= 10:
            processed_records.append(SeqRecord(Seq(re.sub("B|Z|J", "X", str(record.seq))),
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

    # why is this necessary? entries are already 
    # numeric values
    df =df.apply(pd.to_numeric, args=('coerce',))

    n_genomes = len(set(locus2genome.values()))
    n_minimum_genomes = n_genomes-n_missing
    freq_missing = (n_genomes-float(n_missing))/n_genomes
    limit = freq_missing*n_genomes

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


def get_diamond_uniref_top_hits(databases_dir, phylum_filter, n_hits):
    logging.basicConfig(filename='uniref_to_hits.log', level=logging.DEBUG)
    
    logging.info("Parsing fasta...")
    fasta_records = SeqIO.to_dict(SeqIO.parse("diamond_uniref_hit_seqs.fasta", "fasta"))
    logging.info(f'Done. Parsed ({len(fasta_records)}) entries.')
    
    conn = sqlite3.connect("orthology.db")
    cursor = conn.cursor()

    sql2 = 'attach "' + databases_dir + '/ncbi-taxonomy/linear_taxonomy.db" as linear_taxonomy'
    sql4 = 'attach "diamond_uniref.db" as diamond_uniref'

    cursor.execute(sql2,)
    cursor.execute(sql4,)
    conn.commit()
    
    logging.info(f"Retrieve...")
    logging.info(f"Filtering {phylum_filter}...")
    sql_filtered_hits = """SELECT t1.orthogroup, t1.locus_tag, t3.sseqid, t3.hit_count 
       FROM locus_tag2orthogroup t1 INNER JOIN locus_tag2sequence_hash t2 ON t1.locus_tag=t2.locus_tag 
       INNER JOIN diamond_uniref.diamond_uniref t3 ON t2.sequence_hash=t3.qseqid 
       INNER JOIN linear_taxonomy.ncbi_taxonomy t5 ON t3.taxon_id=t5.tax_id 
       WHERE t5.phylum not in ("%s")
       ORDER BY t1.orthogroup, t1.locus_tag, t3.hit_count;"""  % '","'.join(phylum_filter)
    filtered_hits = cursor.execute(sql_filtered_hits,).fetchall()
    logging.info(f"Number of hits: {len(filtered_hits)}")
    logging.info(f"Retrieve top hits excluding phylum of choice")
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
        # keep up to n_hits
        if len(orthogroup2locus2top_hits[orthogroup][locus_tag]) < n_hits:
            orthogroup2locus2top_hits[orthogroup][locus_tag].append(hit_id)

    logging.info(f"Retrieve orthogroup sequences")
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

    logging.info(f"Retrieve top hits sequences from NCBI")
    # for each group, retrieve aa sequence from NCBI
    # write it to fasta file with orthogroup sequences
    # TODO: retrieve from fasta?
    for n,group in enumerate(orthogroup2locus2top_hits):
        if n % 100 == 0:
            print(f"{n} / {len(orthogroup2locus2top_hits)}")
        ortho_sequences = orthogroup2locus_and_sequence[group]
        refseq_sequence_records = []

        # returns a list of list of top hits per locus
        top_hits = [orthogroup2locus2top_hits[group][i] for i in orthogroup2locus2top_hits[group]]
        nr_top_hit = list(set([hit for hits_list in top_hits for hit in hits_list])) # flatten list... yeah, ugly.
        split_lists = chunks(nr_top_hit, 50)

        # multithreading here would be nice to run queries in parallel
        for one_list in split_lists:
            refseq_sequence_records += [fasta_records[acc] for acc in one_list]

        if len(refseq_sequence_records) > 0:
            with open(group + "_nr_hits.faa", 'w') as f:
                for one_ortholog in ortho_sequences:
                    f.write(">%s\n%s\n" % (one_ortholog[0], one_ortholog[1]))
                for record in refseq_sequence_records:
                    name = record.name
                    f.write(">%s\n%s\n" % (name, str(record.seq)))


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


def rename_uniref(uniref_label):
    # UniRef100_UPI0005A8869C 
    return uniref_label.split("UniRef100_")[1]

def setup_diamond_uniref_db(databases_dir, diamond_tsv_files_list):
    logging.basicConfig(filename='uniref_diamond.log', level=logging.DEBUG)
    # load uniref hits into db
    # retrieve taxids from uniref db
    conn = sqlite3.connect("diamond_uniref.db")
    cursor = conn.cursor()

    

    conn_uniref = sqlite3.connect(databases_dir + "/uniref/uniref100.db")
    cursor_uniref = conn_uniref.cursor()


    sql1 = """CREATE TABLE diamond_uniref(hit_count INTEGER, qseqid varchar(200),
        sseqid varchar(200), pident FLOAT, length INTEGER, mismatch INTEGER,
        gapopen INTEGER, qstart INTEGER, qend INTEGER, sstart INTEGER,
        send INTEGER, evalue FLOAT, bitscore FLOAT, taxon_id INTEGER)"""
    cursor.execute(sql1,)
    conn.commit()

    sql = 'insert into diamond_uniref values (?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
    
    diamond_file_list = diamond_tsv_files_list.split(' ')
    for n_file, one_file in enumerate(diamond_file_list):
        diamond_table = pd.read_csv(one_file, sep="\t", header=None)
        
        diamond_table.iloc[:,1] = diamond_table.apply(lambda x: rename_uniref(x[1]), axis=1)
        
        #nr_acc = set(diamond_table.iloc[:,1].to_list())

        accession2taxon_id ={row[0]:row[1] for row in cursor_uniref.fetchall()}
        accession = ''
        count = ''
        # add hit count as first column
        logging.info(f"Iter rows file {n_file}")
        for index, row in diamond_table.iterrows():
            # if new protein, reinitialise the count
            if row[0] != accession:
                accession = row[0]
                count = 1
            else:
                count+=1
            cursor.execute(sql, [count] + row.tolist() + [0])
        conn.commit()

    # index query accession (hash) + hit number
    sql_index_1 = 'create index hitc on diamond_uniref (hit_count);'
    sql_index_2 = 'create index qacc on diamond_uniref (qseqid);'
    sql_index_3 = 'create index sacc on diamond_uniref (sseqid);'
    

    cursor.execute(sql_index_1)
    cursor.execute(sql_index_2)
    cursor.execute(sql_index_3)
    
    conn.commit()
    
    sql = 'select distinct sseqid from diamond_uniref'
    cursor.execute(sql,)
    nr_acc_list = [i[0] for i in cursor.fetchall()]
    nr_acc_set = set(nr_acc_list)
    logging.info(f"Retrieve fasta sequences for {len(nr_acc_list)} entries")
    seqs_uniref = SeqIO.FastaIO.SimpleFastaParser(gzip.open(databases_dir + "/uniref/uniref100.fasta.gz", "rt"))
    fasta_records = []
    acc2taxid = {}
    for n, (title, sequence) in enumerate(seqs_uniref):
        if n % 10000000 == 0:
            logging.info(f'Parsed {n} records -- {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}')
        accession = rename_uniref(title.split(" ")[0])
        if not accession in nr_acc_set:
            continue
        try:
            taxid = re.search("TaxID=(\d+)", title).group(1)
        except:
            taxid = 1
        acc2taxid[accession] = taxid
        fasta_records.append([accession, sequence])
    
    with open("diamond_uniref_hit_seqs.fasta", "w") as outfasta:
        for acc,seq in fasta_records:
            outfasta .write(f">{acc}\n{seq}\n")
    
    logging.info(f"Updating taxids")
    for n, acc in enumerate(acc2taxid):
        if n % 10000 == 0:
            print(f"{n} / {len(acc2taxid)}")
        cursor.execute(f'update diamond_uniref set taxon_id={acc2taxid[acc]} where sseqid="{acc}"')
    sql_index_4 = 'create index txd on diamond_uniref (taxon_id);'
    cursor.execute(sql_index_4)
    conn.commit()
    
    with open("nr_uniref_hits.tab", 'w') as f:
        for acc in nr_acc_list:
            f.write("%s\n" % acc[0])
