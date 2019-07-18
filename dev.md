
## CHLAMDB setup

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

# Web Interface

- [ ] integration of swissprot keywords (possibility to click on it and get complete list of prot, decsription,...)
- [ ] add tcdb classification to "fam" (annotation, phylogenetic profile,...). Include all classification levels (superfamilies,...)
- [ ] idem with EC classification system
- [ ] integrate interpro hierarchy
- [ ] add explanations for hydropathy plot
- [ ] add download page (download KO mapping, interpro mapping, all phylogenies, orthology table,...)
- [ ] search bar: add option to search for TCDB accessions
