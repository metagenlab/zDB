from Bio import Entrez, SeqIO
from Bio.SeqUtils import CheckSum

import sqlite3


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

    fasta_file = "${fasta_file}"

    string_map = open('string_mapping.tab', 'w')
    no_string_mapping = open('no_string_mapping.faa', 'w')

    string_map.write("locus_tag\tstring_id\n")

    records = SeqIO.parse(fasta_file, "fasta")
    no_mapping_string_records = []
    for record in records:
        sql = 'select accession from hash_table where sequence_hash=?'
        cursor.execute(sql, (CheckSum.seguid(record.seq),))
        hits = cursor.fetchall()
        if len(hits) == 0:
            no_mapping_string_records.append(record)
        else:
            for hit in hits:
              string_map.write("%s\t%s\n" % (record.id,
                                              hit[0]))
    SeqIO.write(no_mapping_string_records, no_string_mapping, "fasta")


def convert_gbk_to_faa(gbf_file, edited_gbf):
    print(edited_gbf)
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
