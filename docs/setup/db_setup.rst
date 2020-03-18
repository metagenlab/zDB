PROBLEMS
========

.. warning::

    Latest versions of biopython ignore double quotes in SQL queries.
    Leads to a lot of errors because I use a lot of double quotes to select columns in "where" statments 
    filtering strings. TODO check for a solution. 

    Reason: Rank is a reserved word in MySQL 8.0.2 onwards, but frequently used in the BIoSQL scheme.
    The solution the Biopython community chose is causinge the escape of double quotes in SQL queries.

Quick hack:

Edit file: BioSQL/BioSeqDatabase.py: disable "SET sql_mode='ANSI_QUOTES';"

Edit file: BioSQL/Loader.py: change double quote to "`rank`"

Minimal SETUP
================

.. note::

    This minimal setup include annotations from input GenBank files + comparative data from OrthoFinder.
    It allows to:

    - search and browse annotations
    - perform comparative analysis based on OrthoFinder orthologous groups
    - Plot interactive circular genome maps with circos
    - Plot local alignments relying on orthology informations
    - Blast the database using the various blast flavours

DB setup: download and import biosql schema
----------------------------------------------

.. code-block:: bash

    chlamdb-setup-sqldb.py -d <db_name>

Create a new MySQL database in setup the BioSQL scheme. Assume that mysql password is stored in `SQLWPD` env variable.

Load genbank files
----------------------

.. code-block:: bash

    chlamdb-load-gbk.py -g gbk_edited/*gbk -d 2020_chlamydia_test

- Load gbk sequences & features into the biosql schema
- Create features_* tables for faster access to 
    - locus_tag
    - protein accession
    - gene
    - product
    - seqfeature id (primary key for CDS)

load orthology data
-----------------------

.. code-block:: bash

    chlamdb-load-orthofinder.py -m Orthogroups.txt -d 2020_chlamydia_test

- add “orthogroup” to sqldb *term* table
- add orthogroup for all locus in the seqfeature_qualifier_value table
- create orthogroup prensence absence matrix (comparative_tables_orthology)
	- create locustag2seqfature_id table (custom_tables_locus2seqfeature_id: 	slow, optimization needed)
	- create orthology_detail tablke (deprecoiated but still mandatory)


Setup comparative basic tables
----------------------------------

.. code-block:: bash

	# minimal comparative tables
	# orthology matrix
	chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -o
	
	# orthology matrix: distinguish plasmids from genomes
	chlamdb-setup-comparative-tables.py -d 2020_chlamydia_test -o -a
	identity_closest_homolog

orthogroups consensus annotation
---------------------------------

.. code-block:: bash

    # TODO allow to make statistics for any subset of those data
    chlamdb-get-consensus-orthogroup-annotation.py

Statistics for:

- gene names
- product
- COG
- KO
- domains


Setup old locus table
----------------------

Mandatory by depreciated since synonymous table can be build at the end

Load alignments
-----------------

.. code-block:: bash

    chlamdb-load-alignments.py -a *faa -d 2020_chlamydia_test -c 6

- Calculate identity between pair of sequences
- Create one table/group into orth_<db name>
- Create mean indentity table (obsolete, not working)

TODO: merge individual group tables into one table

6. chlamdb-load-reference-phylogeny.py
--------------------------------------

.. code-block:: bash

    chlamdb-load-reference-phylogeny.py -r core_genome_phylogeny.nwk -d 2020_chlamydia_test -g  ../../data/gbk_edited/*gbk

7. setup taxonomy table
------------------------

.. code-block:: bash

    chlamdb-setup-linear-taxonomy.py -d 2020_chlamydia_test -s linear_taxonomy.db

Might not be strictly necessary (primarily useful to manage the taxnonomy of 
RefSEq and SwissProt hits) but currently necessary for genome statistics.
Bsed on linear_taxonomy.db sqlite database (see snakemake pipeline).

7. chlamdb-setup-genomes-statistics.py
--------------------------------------

.. code-block:: bash

    chlamdb-setup-genomes-statistics.py -d 2020_chlamydia_test

Aptional utilities/annotations
===============================

1. Setup BLAST databases
------------------------



2. Load gene phylogenies
------------------------


.. code-block:: bash

    chlamdb-load-phylogenies.py

3. Load additional annotations
------------------------------

3.1 Load INTERPRO data
+++++++++++++++++++++++

Setup chlamdb-load-hash2locus.py

chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -p # pfam
chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -i # interpro

chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -k # ko


3.2 Load COG data
+++++++++++++++++

chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -c # COG

3.3 Load Kegg data
+++++++++++++++++++



3.4 Load PRIAM data (EC annotation)
+++++++++++++++++++++++++++++++++++

chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -e # EC PRIAM
 

3.5 Load TCDB data (transporters)
+++++++++++++++++++++++++++++++++


3.6 Load psortb data (subcellular localization)
+++++++++++++++++++++++++++++++++++++++++++++++

3.7 Load T3SS effector data
+++++++++++++++++++++++++++


3.8 Load PDB data
++++++++++++++++++


4. Load BLAST results & phylogenies 
------------------------------------

4.1 BLAST vs RefSeq
+++++++++++++++++++

4.2 BLAST vs SwissProt
++++++++++++++++++++++

4.3 Load BBH phylogenies
++++++++++++++++++++++++


5. Add GC content statistics
------------------------------

.. code-block:: bash

	chlamdb-setup-gc-content-tables.py


6. Identification of conserved gene clusters
---------------------------------------------

.. code-block:: bash

	chlamdb-find-conserved-neighborhood.py -d 2019_06_PVC

7. Basic Phylogenetic profiling
--------------------------------

8. add synonymous table (allow to search for RefSeq, Uniprot, uniparc accessions,...)
---------------------------------------------------------------------------------------

- match to uniprot, refseq, accessions to facilitate search



Config optional data
======================

Table with the list of main data. We could add a check that will show an error message is mandatory data is missing.

================================  ================  =============================================
name                              type              status 
================================  ================  =============================================
gbk_files                         mandatory         FALSE
orthology_data                    mandatory         FALSE
orthology_comparative             mandatory         FALSE
orthology_consensus_annotation    mandatory         FALSE
orthogroup_alignments             mandatory         FALSE
old_locus_table                   mandatory         FALSE
reference_phylogeny               mandatory         FALSE
taxnonomy_table                   mandatory         FALSE
genome_statistics                 mandatory         FALSE
BLAST_database                    optional          FALSE
gene phylogenies                  optional          FALSE
interpro_data                     optional          FALSE
interpro_comparative              optional          FALSE
priam_data                        optional          FALSE
priam_comparative                 optional          FALSE
COG_data                          optional          FALSE
COG_comparative                   optional          FALSE
KEGG_data                         optional          FALSE
KEGG_comparative                  optional          FALSE
TCDB_data                         optional          FALSE
psortb_data                       optional          FALSE
T3SS_data                         optional          FALSE
PDB_data                          optional          FALSE
BLAST_refseq                      optional          FALSE 
BLAST_swissprot                   optional          FALSE
BBH_phylogenies                   optional          FALSE
GC_statistics                     optional          FALSE 
gene_clusters                     optional          FALSE 
phylogenetic_profile              optional          FALSE
synonymous_table                  optional          FALSE
================================  ================  =============================================




UNCLEAR PEPENDANCIES
====================

- chlamdb-setup-linear-taxonomy.py


DIVERS & TODO
=============


- Circos plot: possibility to highligh BBH phylum (highlight_BBH= true)
- Taxnonomy circos plots

- If we don’t want to load interpro annotation, add mandatory 	column to orthology_detail 
    - ALTER TABLE orthology_detail ADD TM varchar(10) DEFAULT 'n/a';
    - ALTER TABLE orthology_detail ADD SP varchar(10) DEFAULT 'n/a';

Missing indexes
----------------

- CREATE FULLTEXT INDEX GPF1 ON orthology_detail(gene);
- CREATE FULLTEXT INDEX GPF2 ON orthology_detail(product);
- CREATE FULLTEXT INDEX GPF3 ON orthology_detail(organism);
- CREATE FULLTEXT INDEX GPF4 ON orthology_detail(gene,product,organism);

