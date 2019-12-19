#!/bin/bash -ue
blastp -db /data/databases/uniprot/swissprot/uniprot_sprot.fasta -query chunk_.1 -outfmt 6 -evalue 0.001  -num_threads 4 > chunk_.1.tab
