=======
Methods
=======

---------
Overview
---------

A simplified scheme of ``ChlamDB`` annotation workflow is shown in **Figure 1**. In summary: 

    * annotated genome assemblies integrated into ChlamDB were downloaded from GenBank_ (or 
      RefSeq_ if the GenBank assembly was not annotated).
    * Protein sequences were annotated based on data from muliple source databases and annotation softwares.
    * Protein sequences were clustered into orthologous groups (orthogroups) with OrthoFinder_. 
        * An alignment and a phylogeny was reconstructed for each orthologous group.
    * The closest homologs of each protein were identified in RefSeq_ and SwissProt_
        * A phylogeny including the closest non-PVC homologs was reconstructed for each orthogroup

.. figure:: ../img/workflowV2.svg
    :figclass: align-center

    Figure 1: Simplified annotation workflow.


-----------------
NextFlow pipeline
-----------------

The annotation pipeline summarized in **Figure 1** was fully automated using NextFlow_ . 
The code with detailed software versions and command parameters is available on github here_.

-----------------------------------------------------
Selection of PVC genomes integrated into the database
-----------------------------------------------------

As of Sepetember 2019, they are `2,816` genome assemblies classified as part of 
the PVC superphylum (NCBI Taxid 1783257) on the `NCBI taxonomy website`_. 
``ChlamDB version 2.0`` (June 2019) includes:

    - all complete PVC genomes (June 2019)
    - all draft *Chlamydiae* genomes excluding draft genomes from species of 
      the *Chlamydia* genus with at least 1 complete genome.

.. note::
    Draft genomes of the most studied *Chlamydia* species were excluded because they are highly redundant

.. note::
    Draft genomes of other PVC phlyla were excluded to limit the total number of genomes in the database. 
    Additional genomes of clades of interest might be integrated in future releases of the database.

The list of genomes is available on `ChlamDB home page`_. Complete genomes can be distinguished from draft 
genomes based on the number of contigs ("N. contigs" column) since draft assemblies exhibit more than one contig.


------------------
Protein annotation
------------------


--------------------------------------
Homology search (RefSeq and SwissProt)
--------------------------------------


-----------------------------
Identification of Orthogroups
-----------------------------

Orthogroups were identified with OrthoFinder_. This tools identify orthogroups based on BLASTp (evalue cutoff of 0.001) 
results using the MCL_ software.

.. note::
   ``Orthologs`` are pairs of genes that descended from a single gene in the last common ancestor (LCA) of two species.

.. note::
    An ``orthogroup`` is the group of genes descended from a single gene in the last common ancestor (LCA) of a group of species.
    As gene duplication and loss occur frequently in bacteria, we rarely have exactly one ortholog in each considered genome.

-----------------------------------
Orthogroup alignments & phylogenies
-----------------------------------

---------------------------------------
Calculation of parwise protein identity
---------------------------------------


------------------------------
Circular genome plots (Circos)
------------------------------

------------------
Species phylogeny
------------------


------------------------------------------
Prediction of protein-protein interactions
------------------------------------------


----------------------------------------------------------
Prediction of candidate type III secetion system effectors
----------------------------------------------------------

------------------------------------------
Taxonomic profile of Pfam domains and COGs
------------------------------------------


------------------------------
SQL database and Web Interface
------------------------------




-----------------
Source databases
-----------------


.. table:: Version of source databases used for the annotation
    :width: 800px
    :align: center

    ==================   ======================
    Database name 	     Version
    ==================   ======================
    UniprotKB-UNIPARC    2019.03
    TCDB 	             2019.06
    RefSeq               90
    CDD (COG)            v3.17
    InterProScan         5.35-74.0
    STRING               v11
    PDB                  2019.06
    COG/CDD              v3.17
    KoFam                2019.04.09
    PRIAM                2018.06
    ==================   ======================

-------------------
Software versions
-------------------

.. table:: Version of the main softwares used for protein annotation
    :width: 800 px
    :align: center

    =============   =======
    Software name 	Version
    =============   =======
    FastTree 	    2.1.10
    Diamond      	0.9.24
    OrthoFinder  	2.2.7
    BLAST       	2.7.1
    CheckM      	1.0.12
    KoFamScan    	2019/4/9
    Mafft       	7.407
    =============   =======

-----------------
Code availability
-----------------

=====================================   ===========================================================
Website interface                       https://github.com/metagenlab/chlamdb
Annotation pipeline                     https://github.com/metagenlab/annotation_pipeline_nextflow
Public database download and indexing   https://github.com/metagenlab/annotation_pipeline_snakemake
=====================================   ===========================================================


.. _`NCBI taxonomy website`: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=1783257&lvl=3&p=gcassembly&lin=f&keep=1&srchmode=1&unlock
.. _GenBank: https://www.ncbi.nlm.nih.gov/genbank/
.. _RefSeq: https://www.ncbi.nlm.nih.gov/refseq/
.. _OrthoFinder: https://github.com/davidemms/OrthoFinder
.. _MCL: https://micans.org/mcl/
.. _NextFlow: https://www.nextflow.io/
.. _here: https://github.com/metagenlab/annotation_pipeline_nextflow/blob/master/annotation_pipeline.nf
.. _`ChlamDB home page`: https://chlamdb.ch/#genomes
.. _SwissProt: https://www.uniprot.org/
