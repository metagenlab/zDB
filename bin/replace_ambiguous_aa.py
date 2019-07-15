#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    import re
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    import sys
     
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta', type=str, help="input fasta file")

    args = parser.parse_args()

    records = SeqIO.parse(open(args.input_fasta, "r"), "fasta")
    edited_records = []
    for record in records:
        # replace ambiguous aa by X
        edited_record = SeqRecord(Seq(re.sub("B|Z|J", "X", str(record.seq)), IUPAC.protein),
                                  id=record.id, 
                                  name=record.name,
                                  description=record.description)
        edited_records.append(edited_record)         
    
    SeqIO.write(edited_records, sys.stdout, "fasta")
