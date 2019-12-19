#!/bin/bash -ue
rpsblast -db /data/databases/cdd/Cog -query seq -outfmt 6 -evalue 0.001 -num_threads 4 > blast_result
