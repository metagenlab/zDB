#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import CheckSum

fasta_file = "merged.faa"

nr_fasta = open('nr.faa', 'w')
nr_mapping = open('nr_mapping.tab', 'w')

checksum_nr_list = []

records = SeqIO.parse(fasta_file, "fasta")
updated_records = []

for record in records:

    checksum = CheckSum.crc64(record.seq)
    nr_mapping.write("%s\t%s\n" % (record.id,
                                   checksum))
    if checksum not in checksum_nr_list:
      checksum_nr_list.append(checksum)
      record.id = checksum
      record.name = ""
      updated_records.append(record)

SeqIO.write(updated_records, nr_fasta, "fasta")
