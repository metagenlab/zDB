
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

# TODO

- [ ] extract best swissprot and refseq hits for each sequence (from database), build phylogeny with non-redundant set of BBH (3 refseq + 3 swissprot?)
- [ ] retrieve swissprot annotation score for swissprot hits & uniprotkb exact matches
- [ ] execure gblast vs TCDB for chunks of 1000 sequences (python 2.7)
- [ ] execute T3 effector prediction (3 tools)
- [ ] execute PRIAM for each genome or hmmsearch with PRIAM database?
- [ ] execute macsyfinder for crispr
- [ ] execute macsyfinder for capsule
- [ ] execute macsyfinder for secretion systems
- [ ] retrieve DOORS2 operons
- [ ] run PSORTb version 3.00
- [ ] run genome properties with interprocsn tsv files (https://genome-properties.readthedocs.io/en/latest/index.html, https://www.ebi.ac.uk/interpro/genomeproperties/)
- [ ] http://csbl.bmb.uga.edu/dbCAN/ ==> based on interproscan results

## Ideas

- [ ] KO annotation: compare GHOSTKOALA with eggnog/uniprot annotations
- [ ] execute plast/mmseq2/diamond vs UNIREF90 for chunks of 1000 sequences?
- [ ] metacyc?
- [ ] exact match to COG betfore doing BLAST to reduce computations? => do not blast sequence that are in the COGdb (assign it directly)
- [ ] predict operons when not in doors?

## Notes

- BLAST vs PLAST vs mmseq2 vs diamond?
- install eggnog-mapper?
- execute eggnog-mapper and ghostkoala online?
