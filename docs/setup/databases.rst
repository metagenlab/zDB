Overview
========

This NextFlow_ workflow annotatate multiple bacterial genomes using reference databases, 
infer orthologous groups and build a reference species phylogeny based on core single copy 
orthologs.

See the `documentation of the ChlamDB`_ database for a more detailed description of the workflow.

Input files
===========

- table of genbank/refseq assembly accessions
- ``fna`` files that will be annotated with prokka

Configuration file
==================

.. code-block:: bash

    ////////////////////
    ///// INPUT ////////
    ////////////////////

    // is there any local assemblies in fna that need to be processed
    params.prokka = true
    params.fna_dir = "fna/"

    params.local_sample_sheet = false
    params.ncbi_sample_sheet = "ncbi_assemblies.csv"

    // all databases used by the different scripts should be located in this directory
    params.databases_dir = "/data/databases"



    //////////////////////////////
    ///// PROTEIN ANNOTATION /////
    //////////////////////////////

    params.cog = true
    params.interproscan = true
    params.uniparc = true
    params.tcdb = true
    params.string = true
    params.pdb = true
    params.oma = true
    params.ko = true
    params.tcdb_gblast = false
    params.PRIAM = true
    params.psortb = true
    params.uniprot_goa = true
    params.uniprot_idmapping = true

    params.macsyfinder = true
    params.effector_prediction = true
    params.DeepT3 = true

    // retrieve uniprot annotations & scores 
    params.uniprot_data = false



    //////////////////////////////
    ///// QUALITY CHECK //////////
    //////////////////////////////

    params.checkm = true



    //////////////////////////////
    ///// Orthology //////////////
    //////////////////////////////

    params.orthofinder = true
    params.orthofinder_output_dir = "output"

    // allow (or not) missing data when identifying single copy core genes
    params.core_missing = 0
    params.core_genome_phylogeny_with_fasttree = true

    // build orthogroup phylogeny
    params.orthogroups_phylogeny_with_iqtree = false
    params.orthogroups_phylogeny_with_fasttree = true



    ////////////////////////////////////////////
    ///// HOMOLOGY SEARCH & PHYLOGENIES ////////
    ////////////////////////////////////////////

    // homology search vs RefSeq and SwissProt
    params.blast_swissprot = false
    params.plast_refseq = false
    params.diamond_refseq = true
    params.diamond_refseq_taxonomy = true

    // build phylogenies including closest BBH
    // possibility to filter BBH based on phylum name
    params.refseq_diamond_BBH_phylogeny = true
    params.refseq_diamond_BBH_phylogeny_top_n_hits = 4
    params.refseq_diamond_BBH_phylogeny_phylum_filter = '["Chlamydiae", "Verrucomicrobia", "Planctomycetes", "Kiritimatiellaeota", "Lentisphaerae"]'

    params.interproscan_home = "$params.databases_dir/interproscan/interproscan-latest"



    //////////////////////////////
    ///// CONTAINERS CONFIG //////
    //////////////////////////////

    // Necessary to be able to access the database directory
    // in singularity

    singularity.runOptions = "--bind /data:/data"
    singularity.enabled = true
    singularity.cacheDir = "$baseDir/singularity"

    // The different containers important for the pipeline
    // The container can now be updated by just editing the line here
    // instead of having to do it for every process using the container

    params.chlamdb_container = "metagenlab/chlamdb_annotation:1.0.3"
    params.checkm_container = "metagenlab/checkm:1.0.20"
    params.annotation_container = "metagenlab/annotation-pipeline:1.2"
    params.psort_container = "metagenlab/psort:3.0.6"



    /////////////////////////////
    ///// EXECUTION CONTROL /////
    /////////////////////////////

    process.queue = 'normal'
    process.memory = '2G'
    process.cpus = 40
    params.executor = 'local'

    executor {
    $lsf {
        queueSize = 100
        pollInterval = '30sec'
    }
    $local {
        cpus = 80
        memory = '32 GB'
    }
    }


    conda.cacheDir = "$HOME/miniconda3/nextflow"

    env {
    // necessary to be able to export the python code out
    // of the main nextflow file
    PYTHONPATH = "$baseDir/bin"
    }

Source databases
================

List of public databases and files needed to perform the annotation. 
Most databases are optional depending of the annotation that are configured in the config file.

=============================  =========  ==============  ==============  =================================================
Name                           Mandatory  Size (06.2019)  Size (01.2020)  Path
=============================  =========  ==============  ==============  =================================================
InterproScan                   ?          74G             ?               https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload
UNIPARC                        ?          152G            215G (index?)   ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/
SwissProt/TrEMBL               ?          30G             ?               ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/
UNIPARC interpro annotations   ?          164G            ?               ftp://ftp.ebi.ac.uk/pub/databases/interpro/uniparc_match.tar.gz
Uniprot idmapping file         ?          94G             113G            ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
TCDB                           ?          10M             24M             http://www.tcdb.org/download.php
PDB                            ?          150M            404M            ftp://ftp.wwpdb.org/pub/pdb/derived_data/
RefSeq                         ?          180G            221G (BL+Dmd)   ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/
prot_acc2taxid                 ?          32G             36G             ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
linear NCBI taxonomy           ?          539M            539M            https://github.com/tpillone/ncbitax2lin
kegg KO profiles database      ?          4.8G            5.8G            ftp://ftp.genome.jp/pub/db/kofam/
paperBLAST                     ?          -               2.3G            https://github.com/morgannprice/PaperBLAST#download
STRING                         ?          -               54G+13G         https://stringdb-static.org/download/items_schema.v11.0.sql.gz, https://stringdb-static.org/download/evidence_schema.v11.0.sql.gz
PRIAM                          ?          3.4G            idem            http://priam.prabi.fr/REL_JAN18/Distribution.zip
Uniprot GOA                    no         30G             50G             ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
TOTAL                                     ~765            ?
=============================  =========  ==============  ==============  =================================================


.. note::

    At least ``800Gb`` of disk space needed (TODO: remove temporary/unnecessary files, set them as temp snakemake outputs)


Download & indexing of source databases
======================================

All necessary databases can be downloaded automatically using this `Snakemake workflow`_. 
Since some of these databases are very larges, necessary data are indexed and stored into 
sqlite_ databases. The annotation workflow retrieve annotations either by 

* exact match of amino acid `sequences hashes`_
* using accession numbers
* using various sequence similarity search tools (hmmer, BLAST, Diamond,...)


.. note::
    Please note that downloading and indexing these databases can take between few hours to days.

Snakemake command to download and index most databases except interproscan:

.. code-block:: bash

    snakemake --use-conda --conda-prefix ~/miniconda3/ -j 8 --restart-times 0 
    --snakefile $REPO_PATH/databases_setup/Snakefile_databases 
    update_main


Overview of executed tasks:

.. code-block:: bash

    Job counts:
        count	jobs
        1	download_accession2taxid
        1	download_goa
        1	download_idmapping
        1	download_kofam
        1	download_ncbi_taxonomy
        1	download_paperblast_db
        1	download_paperblast_seqs
        1	download_pdb
        2003	download_refseq_single_file
        1	download_tcdb
        1	download_uniparc
        1	download_uniprotkb_interpro_annotation
        1	extract_ncbi_taxonomy
        1	format_refseq_blast
        1	format_refseq_diamond
        1	formatdb_paperblast
        1	formatdb_pdb
        1	formatdb_tcdb
        1	get_linear_taxonomy
        1	index_accession2taxid
        1	index_id_Source databases
        1	index_interproscan
        1	index_pdb
        1	index_tcdb
        1	index_uniparc
        1	index_uniprotkb_goa
        1	merge_refseq
        1	setup_taxonomy_db
        1	uncompress_koafam
        1	update_main
        1	download_and_uncompress_uniprot_swissprot_fasta
        1	download_swissprot_xml
        1	download_trembl_xml
        1	format_swissprot
        1	index_swissprot_xml
        2038


Snakemake command to download the letest interproscan version:

.. code-block:: bash

    snakemake --use-conda --conda-prefix ~/miniconda3/ -j 10 
    --snakefile $REPO_PATH/databases_setup/Snakefile_databases interproscan

Overview of executed tasks:

.. code-block:: bash

    Job counts:
        count	jobs
        1	copy_missing_executables
        1	download_interproscan
        1	interproscan
        1	uncompress_interproscan
        4

Essential files
+++++++++++++++

====================================================== ==================================
PATH                                                   Description
====================================================== ==================================
interproscan/annotated_uniparc/uniparc_match.db        Index of InterPro annotations for the entire UNIPARC database
uniprot/uniparc/uniparc.db                             Indexed UNIPARC database
uniprot/idmapping/uniprot_sprot_trembl.db              Index UniProt database (annotations)
refseq/merged_refseq.db                                Indexed RefSeq database
ncbi-taxonomy/prot_accession2taxid.db                  Correspondance between RefSeq accession and NCBI taxids
ncbi-taxonomy/linear_taxonomy.db                       Linear NCBI taxnonomy database  
databases/refseq/merged_refseq.dmnd                    RefSeq Diamond database
databases/goa/goa_uniprot.db                           GOA database index
databases/pdb/pdb_seqres.faa.*                         PBD BLASTp database
databases/TCDB/tcdb.faa.*                              TCDB  BLASTp database
databases/kegg/profiles/*                              HMM profiles for KO annotation
====================================================== ==================================


Output files
============

- in order to accelerate protein annotation, identical proteins are annotated only once
- identical protein are identified with their sequence hash through the whole annotation workflow

====================================================== ==================================
PATH                                                   Comment
====================================================== ==================================
data/nr_mapping.tab                                    Correspondance table between locus_tag and sequence-hash
data/nr.faa                                            Non-redundant fasta file (amino acid sequences)
data/filtered_sequences.faa                            Protein sequences larger than 10aa, amibiguous aa replaced by "X"
data/refseq_corresp/refseq_corresp.tab                 Correspondance with RefSeq protein accession & locus_tags
data/checkm/checkm_results.tab                         CheckM results (genome completeness & contamination)
data/gbk_ncbi/*                                        Downloaded GBK files
data/gbk_edited/*                                      Edited gbk files (contigs concatenated, cleaned description, check if all CDS have a locus_tag,... see gbk_check()) 
data/faa_locus/*                                       faa files with CDS locus tag as header
annotation/uniparc_mapping/no_uniparc_mapping.faa      Sequences that could not be mapped to uniparc
annotation/uniparc_mapping/uniparc_mapping.faa         Sequences that could be mapped to uniparc
annotation/uniparc_crossreferences.tab                 Cross reference to many databases extracted from uniparc db                     
====================================================== ==================================

Protein annotation
++++++++++++++++++

- annotations. Generally in tabular format

====================================================== ==================================
PATH                                                   Comment
====================================================== ==================================
annotation/COG/blast_COG.tab                           RPSBLAST results vs cdd COG profiles
annotation/interproscan/*xml,tsv,html/svg              Interproscan results (extracted from db or de novo)
annotation/KO/chunk*                                   KofamScan results 
annotation/paperblast/                                 TODO BLASTp vs DB
annotation/pdb_mapping/blast_results.tab               BLASTp results vs PDB
annotation/PRIAM/sequenceECs.txt                       Ezyme (EC) annotation from PRIAM
annotation/psortb/chunk*                               Subcellular localization from psortB 
annotation/string/string_mapping_PMID.tab              TODO: switch to blast vs db
annotation/T3SS_effectors/BPBAac_results.tab           BPBAac results (T3SS effectors)
annotation/T3SS_effectors/DeepT3_results.tab           DeepT3 results (T3SS effectors)
annotation/T3SS_effectors/effectiveT3_results.tab      effectiveT3 results (T3SS effectors)
annotation/T3SS_effectors/T3_MM_results.tab            T3_MM results (T3SS effectors)
annotation/tcdb_mapping/TCDB_RESULTS*                  TCDB results (transporters)
====================================================== ==================================


Orthology
++++++++++

Orthogroups + pyhlogenies (1/group)
....................................

- fasta files, alignments and phylogenies

====================================================== ==================================
PATH                                                   Comment
====================================================== ==================================
orthology.db                                           Sqlite for fast retrieval and filtering: locus_tag2orthogroup, sequence_hash2aa_sequence, locus_tag2sequence_hash
orthogroups_fasta/OG*faa                               Orthogroup faa files 
orthogroups_alignments/OG*faa                          Orthogroup mafft alignments 
orthogroups_phylogenies_fasttree/*nwk                  Orthogroup FastTree phylogenies 
====================================================== ==================================

Species phylogeny
..................

- core single copy orthogroups, concatenated alignment and species phylogeny

================================================================== ==================================
PATH                                                               Comment
================================================================== ==================================
orthology/core_groups/group*                                       Single copy core groups
orthology/core_alignment_and_phylogeny/msa.faa                     Concatenated alignment of core groups (missing data are replaced by gaps)
orthology/core_alignment_and_phylogeny/core_genome_phylogeny.nwk   Reference species phylogeny reconstructed based on the concatenated alignment with FasTree 
orthology/core_alignment_and_phylogeny/corresp_table.tab           Table of core groups
================================================================== ==================================


Homology search + phylogenies
++++++++++++++++++++++++++++++


SwissProt & RefSeq homologs
...........................

- closest SwissProt & RefSeq homologs, sqlite db for fast data retrieval

========================================================   ==================================
PATH                                                       Comment
========================================================   ==================================
annotation/blast_swissprot/chunk*                          BLASTp results vs SwissProt
annotation/diamond_refseq/chunk*                           Diamond results vs RefSeq
annotation/diamond_refseq/nr_refseq_hits.tab               Non redundant list of hits
annotation/diamond_refseq/diamond_refseq.db                Sqlite database of diamond hits for fast retrival and filtering
annotation/diamond_refseq/refseq_taxonomy.db               Sqlite database with accession to ncbi taxon id correspondance
========================================================   ==================================


Orthogroups + RefSeq BBH pyhlogenies (1/group)
...............................................

- orthogroups fasta with closest RefSeq homologs, alignments and phylogenies

===========================================================   ==================================
PATH                                                          Comment
===========================================================   ==================================
annotation/diamond_refseq_BBH_phylogenies/OG*faa              Faa of each orthologous group + 4 top RefSeq hits of each sequence
orthology/orthogroups_refseq_diamond_BBH_alignments/OG*faa    Mafft alignment of each orthologous group + 4 top RefSeq hits of each sequence
orthology/orthogroups_refseq_diamond_BBH_phylogenies/OG*nwk   FastTree phylogeny of each orthologous groups + 4 top RefSeq hits of each sequence
===========================================================   ==================================


STRING
======

STRING can be downloaded as postgresql database dump. We have to load it to extract:

- PMID associations from textmining (DONE, see ``setup_string_paper_db.py``)
- predicted protein-protein interactions (TODO)
- load into database with ``chlamdb-load-string-pmid-mapping.py``.

PaperBlast
==========

Not yet integrated into the pipeline. Process blast results with ``chlamdb-load-paperblast.py``.

Database setup
==============

Multiple scripts are used to import annotations and compareative data into a MySQL database. `Scripts are here`_
They need to be executed in a specific order (e.g load gbk files and only then the annotations, orthology & phylogenies). 

- TODO: automate setup `with this nextflow workflow`_



.. _Nextflow : https://www.nextflow.io/
.. _`documentation of the ChlamDB` : https://chlamdb.ch/docs/methods/annotation.html
.. _`Snakemake workflow` : https://github.com/metagenlab/databases_setup
.. _sqlite : https://www.sqlite.org/index.html
.. _`sequences hashes` : https://biopython.org/DIST/docs/api/Bio.SeqUtils.CheckSum-module.html
.. _`Scripts are here` : https://github.com/metagenlab/chlamdb/tree/master/db_setup
.. _`with this nextflow workflow` : https://github.com/metagenlab/annotation_pipeline_nextflow/blob/master/chlamdb_setup.nf
