================
Website tutorial
================

.. warning:: 
    This page is under construction. Contact Trestan Pillonel (trestan.pillonel@chuv.ch) if you have any question or suggestion regarding the website or its documentation.

--------------------------------
Search for a protein of interest
--------------------------------

The following protein accessions are recognised:

=============================   =================
Name 	                        Example
=============================   =================
Genbank locus tag 	            wcw_1594_
Genbank protein accession 	    ADI38940.1_
RefSeq locus tag 	            WCW_RS07680_
RefSeq protein accession 	    WP_013182646.1_
UniParc accession 	            UPI0001D5C1DD_
UniProtKB-Accession 	        D6YS95_
UniProtKB-ID 	                D6YS95_WADCW_
=============================   =================

\

\

It's also possible to seach for a gene or product name:

    * mreb_
    * `mreb Waddlia`_
    * `secretion system`_

\

Please note that the search is performed in coding sequence annotations but also in ``EC``, ``Kegg Orthologs``, ``Interpro`` and ``Kegg 
Pathways/Modules`` description fields. The results (if matches were found) are reported in separate tabs (Figure 1).

\

.. figure:: ../img/search.png
    :figclass: align-center
    :width: 90%

    **Figure 1:** Search result for ``secretion system``. Note the presence of multiple tabs with search results in 
    coding sequence annotations ("**locus** tab"), but also in **EC**, **Kegg Orthologs** (KO), 
    **Interpro** and **Kegg Pathways**/**Modules** descriptions.


It's also possible to browse genomes tables from links listed in the `ChlamDB home page`_ (column: Browse online)



-----------------
Species phylogeny
-----------------

----------------------
Orthogroup phylogenies
----------------------

“Putative interactors were predicted in-house from genomic data alone using phylogenetic profiling and investigation of conserved gene neighborhood (see online methods). Phylogenetic profile similarity (pattern of presence/absence of orthologs within the PVC superphylum) was calculated using Euclidean and Jaccard distances. Conserved neighbors were identified by identifying orthologs encoded less than 20 kilobases apart in genomes from different species of the PVC superphylum (Figure 1.H). See (19) and (20) for the rationale justifying use of those two approaches.”


---------------------------------------------------
Search using COG, Pfam, Interpro or KEGG accessions
---------------------------------------------------

Accessions from ``KEGG``, ``COG``, ``Pfam`` and ``InterPro`` can also be searched. 
The result page will report a summary of the entry,  the list of proteins annotated 
with this entry as well as a figure showing the presence/absence of this annotation 
in all genomes included in the database (Figure 2).
\

==================  ==========  =========================================================
Accession type 	    Example 	Description
==================  ==========  =========================================================
KEGG ortholog 	    K00844_ 	hexokinase [EC:2.7.1.1]
COG 	            COG0333_ 	Ribosomal protein L32
PFAM 	            PF06723_ 	MreB/Mbl protein
InterPro            IPR004753_  Cell shape determining protein MreB
KEGG modules 	    M00023_ 	Amino acid metabolism
Kegg pathways 	    map00400_ 	Phenylalanine, tyrosine and tryptophan biosynthesis
==================  ==========  =========================================================

\

\

.. figure:: ../img/K01902_profile.svg
    :figclass: align-center
    :width: 60%

    **Figure 2:** See the `complete profile online`_. 

\

\


.. figure:: ../img/TCA_MAP.svg
    :figclass: align-center
    :width: 90%

------------------------------------------
Taxonomic profile of COGs and Pfam domains
------------------------------------------

----------------
BLAST interface
----------------

A BLAST interface is also available for homology search:

.. figure:: ../img/screenshot_blast.png
    :figclass: align-center

    Figure 1: Simplified annotation workflow.

------------------------
Protein annotation view
------------------------

.. figure:: ../img/locus_page.svg
    :figclass: align-center

    Figure 1: Simplified annotation workflow.

-----------------------------
Orthogroup annotation summary
-----------------------------



---------------------------------------------------------------
Alignments of target genomic regions (from two or more genomes)
---------------------------------------------------------------

.. figure:: ../img/region_align.svg
    :figclass: align-center

    Figure 1: Simplified annotation workflow.

----------------------------------------------------
Whole genomes alignments: interactive circular plots
----------------------------------------------------

.. figure:: ../img/circos_interactive.png
    :figclass: align-center

    Figure 1: Simplified annotation workflow.

-------------------------------------------------------------------
Identification of genes specific to one or more strain(s)/specie(s)
-------------------------------------------------------------------


.. figure:: ../img/extract_orthogroup_page.png
    :figclass: align-center

    Figure 1: Simplified annotation workflow.

--------------------------------------------------------------------------------------------------------
Comparison of the annotation: identification of conserved or clade specific domains/COGs,EC numbers,...
--------------------------------------------------------------------------------------------------------


.. _`ChlamDB home page`: https://chlamdb.ch/#genomes
.. _mreb: https://chlamdb.ch/locusx?accession=mreb
.. _`mreb Waddlia`: https://chlamdb.ch/locusx?accession=mreb+Waddlia
.. _`secretion system`: https://chlamdb.ch/locusx?accession=secretion+system
.. _wcw_1594 : https://chlamdb.ch/locusx?accession=wcw_1594
.. _ADI38940.1 : https://chlamdb.ch/locusx?accession=ADI38940.1
.. _WCW_RS07680 : https://chlamdb.ch/locusx?accession=WCW_RS07680
.. _WP_013182646.1 : https://chlamdb.ch/locusx?accession=WP_013182646.1
.. _UPI0001D5C1DD : https://chlamdb.ch/locusx?accession=UPI0001D5C1DD
.. _D6YS95 : https://chlamdb.ch/locusx?accession=D6YS95
.. _D6YS95_WADCW : https://chlamdb.ch/locusx?accession=D6YS95_WADCW
.. _K00844 : https://chlamdb.ch/locusx?accession=K00844
.. _COG0333 : https://chlamdb.ch/locusx?accession=COG0333
.. _PF06723 : https://chlamdb.ch/locusx?accession=PF06723
.. _IPR004753 : https://chlamdb.ch/locusx?accession=IPR004753
.. _M00023 : https://chlamdb.ch/locusx?accession=M00023
.. _map00400 : https://chlamdb.ch/locusx?accession=map00400
.. _`complete profile online` : https://chlamdb.ch/locusx?accession=K01902#tab3

