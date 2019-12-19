#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum

conn = sqlite3.connect("/data/databases/uniprot/uniparc/uniparc.db")
cursor = conn.cursor()

fasta_file = "merged_nr.faa"

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
    sql = 'select t1.uniparc_id,uniparc_accession,accession,taxon_id,description, db_name, status from uniparc_accession t1 inner join uniparc_cross_references t2 on t1.uniparc_id=t2.uniparc_id inner join crossref_databases t3 on t2.db_id=t3.db_id where sequence_hash=?'
    cursor.execute(sql, (CheckSum.seguid(record.seq),))
    hits = cursor.fetchall()
    if len(hits) == 0:
        no_mapping_uniparc_records.append(record)
        no_mapping_uniprot_records.append(record)
    else:
        mapping_uniparc_records.append(record)
        all_status = [i[6] for i in hits]
        if 1 in all_status:
            status = 'active'
        else:
            status = 'dead'
        uniparc_map.write("%s\t%s\t%s\t%s\n" % (record.id,
                                               hits[0][0],
                                               hits[0][1],
                                               status))
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
