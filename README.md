
# ChlamDB django app


```

update seqfeature_qualifier_value set value="RB11436b" where seqfeature_id=147658 and value="RB11436";

OK chlamdb-load-gbk.py -g *gbk -d 2019_06_PVC
OK chlamdb-load-orthofinder.py -m Orthogroups.txt -d 2019_06_PVC
chlamdb-setup-old_locus-table.py -d 2019_06_PVC
chlamdb-setup-genomes-statistics.py -d 2019_06_PVC

chlamdb-load-alignments.py -a *faa -d 2019_06_PVC
chlamdb-load-reference-phylogeny.py -r core_genome_phylogeny.nwk -d 2019_06_PVC -g ../../data/gbk_edited/*gbk
chlamdb-load-swissprot-homology-search.py -i chunk_.*.tab -d 2019_06_PVC -t -p 2 -l -u ../../data/nr_mapping.tab


chlamdb-load-COG.py -i blast_COG.tab -d 2019_06_PVC -u ../../data/nr_mapping.tab -cc cog_corresp.tab -cl cog_length.tab
chlamdb-load-interproscan.py -i *tsv -d 2019_06_PVC -u ../../data/nr_mapping.tab
chlamdb-load-KO.py -k chunk*.tab -d 2019_06_PVC -c ../../data/nr_mapping.tab
chlamdb-load-PRIAM.py -i sequenceECs.txt -d 2019_06_PVC -c ../../data/nr_mapping.tab


chlamdb-load-uniprot-annotations.py -d 2019_06_PVC -um uniprot_mapping.tab -d uniprot_data.tab -hm ../../data/nr_mapping.tab

chlamdb-setup-linear-taxonomy.py -d 2019_06_PVC -s linear_taxonomy.db
chlamdb-setup-gc-content-tables.py -d 2019_06_PVC


```
