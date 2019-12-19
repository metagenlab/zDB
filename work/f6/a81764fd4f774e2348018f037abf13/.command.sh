#!/usr/bin/env python


from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum

conn = sqlite3.connect("/data/databases/oma/oma.db")
cursor = conn.cursor()

fasta_file = "merged_nr.faa"

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
    else:
        for hit in hits:
          oma_map.write("%s\t%s\n" % (record.id,
                                              hit[0]))


SeqIO.write(no_oma_mapping_records, no_oma_mapping, "fasta")
