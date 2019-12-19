#!/bin/bash -ue
diamond blastp -p 4 -d /data/databases/refseq/merged_refseq.dmnd -q chunk_.1 -o chunk_.1.tab --max-target-seqs 200 -e 0.01 --max-hsps 1
