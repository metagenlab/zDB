
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
- execute interproscan
- mapping to string, oma, pdb, tcdb databases (exact matches)
- get PMID mapping from string

# refseq reference/representative genomes

- download from RefSeq database
- annotate each genome with InterproScan: TODO: use exact match to uniparc to get precomputed annotations, only execute interproscan for unannotated proteins

# DONE


- [X] parallelized BLASTp for orthofinder
- [X] orthofinder
- [X] alignments of each orthogroup with mafft
- [X] reconstruction of a phylogeny for orthogroups with more than 3 sequences with iq-tree
- [X] identify single copy orthogroups (accept missing data)
- [X] write core orthogroups with accessions as headers
- [X] concatenate core orthogroups
- [X] build phylogeny with fasttree
- [X] remove reducancy before annotation (use sequence hash as seqid)
- [X] execute blast vs SwissProt for chunks of 1000 sequences
- [X] execute rpsblast vs COG for chunks of 1000 sequences
- [X] get mapping to uniparc
- [X] get exact match to eggnog/string /home/tpillone/work/dev/bio_databases/tri/sqlite/seq_db_eg.db
- [X] get exact match to oma
- [X] get exact match to PDB
- [X] get PMID associations from string database
- [X] execute plast RefSeq for chunks of 1000 sequences
- [X] execute diamond RefSeq for chunks of 1000 sequences
- [X] execute interproscan by chunks of 1000 sequences
- [X] execure gblast vs TCDB for chunks of 1000 sequences (only unmapped sequences): todo: solve dependancies/env variables

# TODO

- [ ] extract best swissprot and refseq hits for each sequence (from database), build phylogeny with non-redundant set of BBH (3 refseq + 3 swissprot?)
- [ ] retrieve GO annotations from uniprotKB GOA (exact match or best diamond/plast hit if no exact match)
- [ ] retrieve uniprot annotation from uniprotKB (exact match or best diamond/plast hit if no exact match)
- [ ] retrieve uniprot annotation score for matched uniprotkb entries

- [ ] execute T3 effector prediction: BPBAac
- [ ] execute T3 effector prediction: T3_MM
- [ ] execute T3 effector prediction: effective

- [ ] execute macsyfinder for crispr
- [ ] execute macsyfinder for capsule
- [ ] execute macsyfinder for secretion systems

- [ ] execute PRIAM for each genome or hmmsearch with PRIAM database?

- [ ] retrieve DOORS2 operons
- [ ] predict operons when not in doors? operon-mapper?

- [ ] run genome properties with interprocsn tsv files (https://genome-properties.readthedocs.io/en/latest/index.html, https://www.ebi.ac.uk/interpro/genomeproperties/)
- [ ] http://csbl.bmb.uga.edu/dbCAN/ ==> based on interproscan results

- [ ] PSORTb version 3.00

- [ ] get kegg annotation from eggnog, unpiprotKB extact or best match, eggnog

## Ideas

- [ ] IS annotation with https://github.com/xiezhq/ISEScan (https://github.com/emrobe/SfB-course/blob/master/Deployment_wrapper.sh) or https://github.com/thanhleviet/ISfinder-sequences
- [ ] KO annotation: compare GHOSTKOALA with eggnog/uniprot annotations
- [ ] execute plast/mmseq2/diamond vs UNIREF90 for chunks of 1000 sequences?
- [ ] metacyc?
- [ ] exact match to COG betfore doing BLAST to reduce computations? => do not blast sequence that are in the COGdb (assign it directly)
- [ ] predict operons when not in doors?

## Notes

- ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/bacsu.txt
- ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/dbxref.txt

- BLAST vs PLAST vs mmseq2 vs diamond: statisics and best hits comparisons
- install eggnog-mapper?
- execute eggnog-mapper and ghostkoala online?
