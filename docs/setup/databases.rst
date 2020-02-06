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

=============================  =========  =================================================
Name                           Mandatory  Path
=============================  =========  =================================================
interproscan                   ?          https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload
SwissProt/TrEMBL               ?          ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/
Uniprot idmapping file         ?          ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
UNIPARC                        ?          ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/
TCDB                           ?          http://www.tcdb.org/download.php
PDB                            ?          ftp://ftp.wwpdb.org/pub/pdb/derived_data/
RefSeq                         ?          ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/
prot_acc2taxid                 ?          ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
linear NCBI taxonomy           ?          https://github.com/tpillone/ncbitax2lin
kegg KO profiles database      ?          ftp://ftp.genome.jp/pub/db/kofam/
UNIPARC interpro annotations   ?          ftp://ftp.ebi.ac.uk/pub/databases/interpro/uniparc_match.tar.gz
paperBLAST                     ?          https://github.com/morgannprice/PaperBLAST#download
STRING                         ?          https://stringdb-static.org/download/items_schema.v11.0.sql.gz, https://stringdb-static.org/download/evidence_schema.v11.0.sql.gz
PRIAM                          ?          http://priam.prabi.fr/REL_JAN18/Distribution.zip
=============================  =========  =================================================



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


Output files
============






.. _Nextflow : https://www.nextflow.io/
.. _`documentation of the ChlamDB` : https://chlamdb.ch/docs/methods/annotation.html
.. _`Snakemake workflow` : https://github.com/metagenlab/databases_setup
.. _sqlite : https://www.sqlite.org/index.html
.. _`sequences hashes` : https://biopython.org/DIST/docs/api/Bio.SeqUtils.CheckSum-module.html
