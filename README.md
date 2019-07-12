
# annotation pipeline (chlamdb-type database)

- automatic download from GenBank (input: assembly table)
- check gbk format for biosqldb setup (concatenate draft genomes if necessary)
- blast each pairs of proteomes
- run orthoFinder
- extract orthogroups fasta
- align each orthogroup with mafft
- build phylogenies with iq-tree or fasttree
- identify core single copy orthogroups
- concatenate core orthogroup alignments & build a reference phylogeny with FastTree
- rps-BLAST to COG database (cdd PSSM)
- BLASTp to swissprot database
- diamond or plast to RefSEq (top 100 hits)
- annotate with interproscan
- execute gBLAST for TCDB annotation (transporters)
- execute KofamScan to get KO annotation
- mapping to string, oma, pdb, tcdb databases (exact matches)
- get PMID mapping from string


# DONE

- [X] input samples table1: genbank assembly accession
- [X] input samples table2: path to local gbk files
- [X] gbk_check (add unique taxid for each genome, concatenate contigs)
- [X] parallelized BLASTp for orthofinder
- [X] orthofinder
- [X] alignments of each orthogroup with mafft
- [X] reconstruction of a phylogeny for orthogroups with more than 3 sequences with iq-tree
- [X] identify single copy orthogroups. Accept missing data (as parameter, default=0 missing data)
- [X] write core orthogroups with accessions as headers
- [X] concatenate core orthogroups
- [X] build phylogeny with fasttree
- [X] remove aa sequences reducancy before annotation (use sequence hash as seqid)
- [X] execute blast vs SwissProt for chunks of 1000 sequences
- [X] execute rpsblast vs COG for chunks of 1000 sequences
- [X] execute plast RefSeq for chunks of 1000 sequences
- [X] execute diamond RefSeq for chunks of 1000 sequences
- [X] execute interproscan by chunks of 1000 sequences
- [X] execute kofamscan by chunks of 1000 sequence (https://www.genome.jp/tools/kofamkoala/)
- [X] get mapping to uniparc (exact match)
- [X] get exact match to eggnog/string
- [X] get exact match to oma
- [X] get exact match to PDB
- [X] get PMID associations from string database for mapped seqs
- [X] execute gblast vs TCDB for chunks of 1000 sequences (only for unmapped sequences): todo: solve dependancies/env variables
- [X] interproscan: analyse separately uniparc matches from non matched entries (allocate more ressources to the later): no match: more ram (16go), and cores (8), no lookup service
- [X] extract best refseq (and swissprot?) hits for each sequence (from database), build phylogeny with non-redundant set of BBH (x top seqs refseq + y top seqs swissprot as parameters)
  - [X] get taxid of hits => make protein_accession2phylum sqlite db
  - [X] setup db (single process)
    - [X] load blast results into sqlite db
    - [X] load locus tag 2 hash table
    - [X] load orthology data into sqlite db (locus2orthogroup)
  - [X] keep X hits/sequence filter out Chlamydiae/Planctomycetes,... hits (attach linear taxonomy)
  - [X] retrieve refseq hits aa sequences
  - [X] generate one fasta files/orthogroup
  - [X] align with mafft
  - [X] build phylogeny
- [X] refseq taxonomy: use downloaded accession2taxonomy + indexed refseq fasta instead of web queries
- [X] retrieve uniprot annotation score for matched uniprotkb entries
- [X] execute checkM

# TODO

## CHLAMDB setup

## priority 1

### setup

- [X] setup biosqldb schema
- [ ] setup kegg tables (ENZYME database)
  - [ ] separate script into multiple scripts and deal with incomplete tables
  - [X] move to separate script (not in load script anymore)
  - [X] ko2annotation table
    - [ ] redundancy with ko2module and ko2pathway, can be simplified and accelerated
  - [X] setup enzyme_dat and enzyme tables
  - [X] setup pathway table
  - [X] setup ko2pathway
  - [X] setup module table
  - [X] setup ko2module
- [X] setup COG tables (COG database)
- [X] setup linear taxonomy (with taxid for each rank see virulence db setup)
- [ ] setup TCDB tables
  - [X] download faa: setup_transporter table, tc_table and uniprot_table (TODO: rename table)
  - [ ] update table organization (dev), not all entries have uniprot accessions

### annotation results

- [X] load genomes
  - [X] setup features table
  - [X] setup seqfeature_id2locus_table
- [X] load orthofinder results  
- [X] load interproscan results
  - [X] add TM and SP columns to orthology_detail legacy table
- [X] load orthogroup alignments
- [X] get identity closest homolog table
- [X] get average identity table
- [X] setup genome statistics table
- [X] load COG hits
- [X] load KO hits
  - [X] load legacy tables
- [X] setup pairwise BBH identity tables
- [X] lead orthogroup alignments(identity matrices)
- [X] load orthogroup phylogenies
- [X] load orthogroup BBH phylogenies
- [X] load reference phylogeny
- [X] get genome table (homepage)
- [X] get conserved neighborhood
- [X] load uniprot annotations
- [X] load blast swissprot results
  - [X] download taxonomy-description information(s)
- [X] load blast refseq results
    - [X] load refseq taxonomy table
- [X] load TCDB annotations
- [X] get phylo profile
  - [ ] setup core_orthogroups_identity_msa_
- [X] legacy COG table
- [ ] legacy locus2EC table
- [X] legacy PFAM table
- [X] setup blast databases
- [X] phylogenetic profiles
- [X] use celery for circos_main view
- [X] load PMID string mapping
- [X] load checkM
- [ ] load T3SS effector predictions
- [ ] get effectiveT3 "eukaryote" domains
- [ ] load pdb best hits
- [ ] load cross-references
- [ ] update PMID
- [ ] deal with search for KO, KEGG, IP absent from genomes included in the database
- [ ] show confidence scores for PDB, KEGG, COG,... (identity, score, evalue,...)
- [ ] add uniprot proteome column (+ percent overlap)
- [ ] load RBBH data (for identity distribution plots)
- [ ] get species table ==> display on homepage
- [ ] check indexes
- [ ] update browse genome view
- [ ] add uniprot kewords to locus page

## priority 2

- [ ] switch to new COG tables
- [ ] switch to new KEGG tables
- [ ] setup NOG_table_v5 and NOG_members_v5
- [ ] setup interpro master table
- [ ] setup pfam master table
- [ ] load DOOR2 data
    - [ ] accession table
    - [ ] operons data

## Annotation pipeline

## priority 1

- [X] get proteome match
- [ ] retrieve GO annotations from uniprotKB GOA (exact match or best diamond/plast hit if no exact match?)
- [X] execute T3 effector prediction: BPBAac
- [X] execute T3 effector prediction: T3_MM
- [X] execute T3 effector prediction: effectiveT3
- [X] execute T3 effector prediction: DeepT3
- [X] BLASTp vs pdb
- [ ] get cross-references from uniprot IdMapping? Or from uniprot db itself? (cross references from indexed uniprotKB xml)
  - [ ] priority to uniprot entries from corresponding proteome (otherwise based on exact match)
  - [ ] multiple match case? get entire proteomes to get the correct mapping between locus_tag and un iprot entries: https://www.ebi.ac.uk/proteins/api/doc/#!/uniparc/getByProteomeId
- [ ] get Refseq protein ID and locus_tag and match to new locus tags in index

## priority 2

- [ ] get STRING associations if available

- [ ] Inc prediction based on bi-lobbed hydrophobic domains

- [ ] VF annotation with all databases

- [ ] retrieve DOORS2 operons

- [ ] predict operons with cluster_finder.pl and operon_finder.pl

- [ ] mapping to uniprot proteomes: https://www.ebi.ac.uk/proteins/api/doc/#!/proteomes/search
- [ ] curl -X GET --header 'Accept:application/json' 'https://www.ebi.ac.uk/proteins/api/proteomes?offset=0&size=100&xref=GCA_000068525.2'

- [ ] match to uniparc: use https://www.ebi.ac.uk/proteins/api/doc/#!/uniparc/getBySequence
- [ ] get entire proteomes: https://www.ebi.ac.uk/proteins/api/doc/#!/uniparc/getByProteomeId

## development

- [ ] integration of swissprot keywords (possibility to click on it and get complete list of prot, decsription,...)

- [ ] (retrieve uniprot annotation from uniprotKB (exact match or best diamond/plast hit if no exact match?))

- [ ] Resfams annotation

- [ ] execute macsyfinder for crispr
- [ ] execute macsyfinder for capsular genes
- [ ] execute macsyfinder for secretion systems

- [ ] execute PRIAM for each genome or hmmsearch with PRIAM database?

- [ ] predict operons when not in doors? operon-mapper?

- [ ] run genome properties with interprocsn tsv files (https://genome-properties.readthedocs.io/en/latest/index.html, https://www.ebi.ac.uk/interpro/genomeproperties/)

- [ ] Carbohydrate-Active Enzyme database: http://csbl.bmb.uga.edu/dbCAN/

- [ ] PSORTb version 3.00

- [ ] get kegg annotation from eggnog, unpiprotKB extact or best match, eggnog

- [ ]  PMID mapping with ppaperBLAST http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi

# Web Interface

- [ ] add tcdb classification to "fam" (annotation, phylogenetic profile,...). Include all classification levels (superfamilies,...)
- [ ] idem with EC classification system
- [ ] integrate interpro hierarchy
- [ ] add explanations for hydropathy plot
- [ ] add download page (download KO mapping, interpro mapping, all phylogenies, orthology table,...)
- [ ] search bar: add option to search for TCDB accessions

## Ideas

- [ ] check if BioVx tools could be used to identitify candidate Incs
- [ ] https://www.uniprot.org/help/api_idmapping
- [ ] IS annotation with https://github.com/xiezhq/ISEScan (https://github.com/emrobe/SfB-course/blob/master/Deployment_wrapper.sh) or https://github.com/thanhleviet/ISfinder-sequences
- [ ] KO annotation: compare GHOSTKOALA with eggnog/uniprot annotations
- [ ] execute plast/mmseq2/diamond vs UNIREF90 for chunks of 1000 sequences?
- [ ] metacyc?
- [ ] exact match to COG betfore doing BLAST to reduce computations? => do not blast sequence that are in the COGdb (assign it directly)
- [ ] predict operons when not in doors?
- [ ] orthoDB
- [ ] ICEberg
- [ ] viruses: vFam (http://derisilab.ucsf.edu/software/vFam/)
- [ ] VirSorter (http://merenlab.org/2018/02/08/importing-virsorter-annotations/)
- [ ] https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/09-KEGG_programming.html
- [ ] read BioV (gblat3) manual: transporter annotation, programmatic Interface with TCDB,...
- [ ] transporters https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0192851
- [ ] transporters https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0208151, https://www.nature.com/articles/nbt.4024
- [ ] use GTDB instead of NCBI taxonomy?
- [ ] phylogenetic profiling for pfam/interpro profiles with RefSeq data
- [ ] MCL BLAST graph vizualisation
- [ ] standardize annotation within orthogroups: if sufficient reciprocal seq coverage (BLAST align), expand KO, COG and EC annotations   
- [ ] curated blast db: https://msystems.asm.org/content/4/2/e00072-19

## Notes

- ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/bacsu.txt
- ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/dbxref.txt

- BLAST vs PLAST vs mmseq2 vs diamond: statisics and best hits comparisons
- install eggnog-mapper?
- HHBLITS

# refseq reference/representative genomes

- [X] download representative/reference genomes from RefSeq database
- [X] annotate each genome with InterproScan: use exact match to uniparc to get precomputed annotations
- [X] annotate representative genomes with kofamscan => phylogenetic profile of KO, modules, pathways,...
- [X] execute local interproscan for unannotated proteins
- [X] execute rpsblast COG cdd 
- alternative option: work with uniprotkb proteomes 1) exclude anomalous proteomes based on refseq data 2) get species taxid for each proteome 3) remove redudancy (keep one prepresentative per species taxid) 4) retrieve interpro annotation from interproscan uniparc annotations
- [ ]make stats from ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt: superkingdom, annotated genomes,...
- [ ] use GTDB rather than NCBI taxonomy? ==> does not include eukaryotes
  - [ ] see http://annotree.uwaterloo.ca/app/#/?qtype=pfam&qstring=PF00617&eval=0.00001
  - [ ] get leaf2phylum & leaf2order & leaf2class
  - [ ] extract subset of their reference phylogeny
  - [ ] get phylum tree & order tree & class tree (keep only lineages containing representative/reference genomes) (see how many different phylum, class, orders)
