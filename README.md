
# ChlamDB django app


```

update seqfeature_qualifier_value set value="RB11436b" where seqfeature_id=147658 and value="RB11436";

OK chlamdb-load-gbk.py -g *gbk -d 2019_06_PVC
OK chlamdb-load-orthofinder.py -m Orthogroups.txt -d 2019_06_PVC
OK chlamdb-setup-old_locus-table.py -d 2019_06_PVC
OK chlamdb-setup-genomes-statistics.py -d 2019_06_PVC
OK chlamdb-load-hash2locus.py -u nr_mapping.tab -d 2019_06_PVC


OK chlamdb-load-reference-phylogeny.py -r core_genome_phylogeny.nwk -d 2019_06_PVC -g ../../data/gbk_edited/*gbk

# load interproscan results
OK chlamdb-load-interproscan.py -i *tsv -d 2019_06_PVC -u ../../data/nr_mapping.tab

# add legacy table
OK chlamdb-load-interproscan.py -i *tsv -d 2019_06_PVC -u ../../data/nr_mapping.tab -l

# add add_SP_TM to orthology_ table
OK chlamdb-load-interproscan.py -a -d 2019_06_PVC

# load COG and legacy table
OK chlamdb-load-COG.py -i blast_COG.tab -d 2019_06_PVC -u ../../data/nr_mapping.tab -cc cog_corresp.tab -cl cog_length.tab -l

OK chlamdb-load-KO.py -k chunk*.tab -d 2019_06_PVC -c ../../data/nr_mapping.tab

chlamdb-load-uniprot-annotations.py -d 2019_06_PVC -um uniprot_mapping.tab -d uniprot_data.tab -hm ../../data/nr_mapping.tab


OK chlamdb-load-alignments.py -a *faa -d 2019_06_PVC -c 100
chlamdb-load-swissprot-homology-search.py -i chunk_.*.tab -d 2019_06_PVC -t -p 2 -l -u ../../data/nr_mapping.tab

chlamdb-load-PRIAM.py -i sequenceECs.txt -d 2019_06_PVC -c ../../data/nr_mapping.tab

# comparative tables
chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -o # orthogroup
chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -c # COG
chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -p # pfam
chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -i # interpro
chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -k # ko

chlamdb-setup-linear-taxonomy.py -d 2019_06_PVC -s linear_taxonomy.db
chlamdb-setup-gc-content-tables.py -d 2019_06_PVC


```
