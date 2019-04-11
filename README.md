
# annotation pipeline (chlamdb-type database)

- blast each pairs of proteomes
- run orthoFinder
- extract orthogroups fasta
- align each orthogroup with mafft
- build phylogenies with iq-tree
- identify core single copy orthogroups
- concatenate core orthogroup alignments & build a reference phylogeny with FastTree
- rps-BLAST to COG database (cdd PSSM)
- BLASTp to swissprot database
- execute interproscan TODO: use exact match to uniparc to get precomputed annotations, only execute interproscan for unannotated proteins
- mapping to string, oma, pdb, tcdb databases (exact matches)
- PDB id mapping from string

# refseq reference/representative genomes

- download from RefSeq database
- annotate each genome with InterproScan: TODO: use exact match to uniparc to get precomputed annotations, only execute interproscan for unannotated proteins

# DONE

- [X] orthofinder
- [X] alignments of each orthogroup with mafft
- [X] reconstruction of a phylogeny for orthogroups with more than 3 sequences with iq-tree
- [X] identify single copy orthogroups (accept missing data)
- [X] write core orthogroups with accessions as headers
- [X] concatenate core orthogroups
- [X] build phylogeny with fasttree
- [X] execute interproscan for each faa
- [X] execute blast vs SwissProt for chunks of 1000 sequences
- [X] execute blast vs COG for chunks of 1000 sequences
- [X] execute profile COG search with CDD database
- [X] get mapping to uniparc
- [X] exact match to eggnog/string /home/tpillone/work/dev/bio_databases/tri/sqlite/seq_db_eg.db
- [X] exact match to oma
- [X] exact match to PDB
- [X] get PMID associations from string database

# TODO

- [ ] exact match to COG betfore doing BLAST to reduce computations? => do not blast sequence that are in the COGdb (assign it directly)

- [ ] interproscan annotation 1) from uniparc annotations 2) from local interproscan analysis for unmapped sequences

- [ ] retrieve swissprot annotation score

- [ ] get DOORS operons
- [ ] predict operons when not in doors? operon-mapper?

- [ ] (execute plast/mmseq2/diamond vs RefSeq for chunks of 1000 sequences)
- [ ] extract best swissprot and refseq hits for each sequence (from database), build phylogeny with non-redundant set of BBH (3 refseq + 3 swissprot?)
- [ ] execute plast/mmseq2/diamond vs UNIREF90 for chunks of 1000 sequences?
- [ ] execure gblast vs TCDB for chunks of 1000 sequences (python 2.7)
- [ ] execute T3 effector prediction (3 tools)
- [ ] execute PRIAM for each genome or hmmsearch with PRIAM database?
- [ ] execute macsyfinder for crispr
- [ ] execute macsyfinder for capsule
- [ ] execute macsyfinder for secretion systems

- [ ] integrate genome properties (https://genome-properties.readthedocs.io/en/latest/index.html, https://www.ebi.ac.uk/interpro/genomeproperties/)
- [ ] http://csbl.bmb.uga.edu/dbCAN/ ==> based on interproscan results
- [ ] metacyc?
- [ ] KO annotation: compare GHOSTKOALA with eggnog/uniprot annotations

* Notes

- BLAST vs PLAST vs mmseq2 vs diamond?
- execute eggnog-mapper and ghostkoala online?
