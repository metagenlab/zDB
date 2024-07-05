=======
Methods
=======

---------
Overview
---------

A simplified scheme of ``zDB`` annotation workflow is shown in **Figure 1**. In summary: 

    * Input consists of annotated genome assemblies. These could be downloaded from GenBank_ or collected and assembled by the user.
    * Protein sequences are annotated based on data from multiple source databases and annotation softwares.
    * Protein sequences are clustered into orthologous groups (orthogroups) with OrthoFinder_.
        * An alignment and a phylogeny is reconstructed for each orthologous group.
    * The closest homologs of each protein are identified in RefSeq_ and SwissProt_
    * Data is then integrated into a custom SQL database derived from the bioSQL_ relational model.


.. figure:: img/workflowV2.png
    :figclass: align-center

    **Figure 1**: Simplified annotation workflow.


The pipelines to setup the reference databases and perform the analysis outlibed in **Figure 1** are fully automated using NextFlow_ (except for the setup of the RefSeq_ database) and can be run simply with the `zdb setup` and `zdb run` commands.


------------------
Protein annotation
------------------

Protein sequences were annotated using data from multiple databases and various bioinformatic softwares (*Figure 1*) as detailed below:
    * ``COG`` annotation is done with rpsblast_ (with an e-value cutoff of *0.001*) using `position-specific scoring matrix`_ (PSSM) from the `Conserved Domain Database`_ (CDD)
    * Proteins are assigned to `Kegg Orthologs`_ (KO) numbers with KofamScan_ (default parameters) and KOfam_ (a database of Hidden Markov Models [HMM] of `KEGG Orthologs`_). KofamScan relies on HMMER_ to search the HMM database.
        * The `KEGG API`_ is used to retrieve the annotation of ``Kegg Orthologs`` and mapping to ``KEGG modules``, ``KEGG pathways``.

    * `PfamScan.pl` is used to scan the set of non-redundant protein sequences against the Pfam_ database to add protein family annotations.
    * Antimicrobial resistance genes are prediced using the AMRFinderPlus_ software
    * Virulence factors are predicted with a protein BLAST_ of the non-redundant protein sequences against the `Virulence Factor Database`_.


-----------------------------------------
_`Homology search` (RefSeq and SwissProt)
-----------------------------------------

* The closest identifiable RefSeq_ homologs are searched with Diamond_ (parameters: ``--max-target-seqs`` 200 ``-e`` 0.01 ``--max-hsps`` 1). The 200 first hits were retained for each protein.

* The closest identifiable SwissProt_ homologs were searched with protein BLAST_ (parameters: ``-evalue`` 0.001).


--------------------------------
_`Identification of Orthogroups`
--------------------------------

Orthogroups (or orthologous groups) were identified with OrthoFinder_. This tools identifies orthogroups based on protein BLAST_ (parameters: ``-evalue`` 0.001) results using the MCL_ clustering software.

.. note::
   ``Orthologs`` are pairs of genes that descended from a single gene in the last common ancestor (LCA) of two species.

.. note::
    An ``orthogroup`` is the group of genes descended from a single gene in the last common ancestor (LCA) of a group of species.
    As gene duplication and loss occur frequently in bacteria, we rarely have exactly one ortholog in each considered genome.


--------------------------------------------------------
Orthogroup multiple sequence alignments & phylogenies
--------------------------------------------------------

Protein sequences of each orthogroup are aligned with MAFFT_ (default parameters). A phylogeny is then reconstructed for each orthogroup of three or more sequences. The phylogeny is reconstructed with FastTree_ with default parameters. The node support values are not traditionnal `boostrap support values`_. FastTree_ uses the Shimodaira-Hasegawa test with 1,000 bootstrap replicates to quickly estimate the reliability of each split in the tree. Values higher than **0.95** can be considered as "strongly supported".


---------------------------------------
Phylogeny including top RefSeq hits
---------------------------------------

A second phylogeny is reconstructed for each orthogroup. This phylogeny includes the 4 best non-PVC RefSeq hits of each protein (see the `Homology search`_ paragraph).
First, the NCBI taxon ID of each RefSeq hit is retrieved using the ``prot.accession2taxid`` mapping file available from the `NCBI taxonomy ftp website`_. PVC hits were removed and the amino acid sequence of the 4 best non-PVC hits is retrieved from the NCBI using the `biopython interface to Entrez`_. Protein sequences of each orthogroup + RefSeq homologs are aligned with MAFFT_ (default parameters) and the phylogeny is reconstructed with FastTree_ (default parameters).

\

This phylogeny allows users to check whether PVC proteins form a monophyletic group or if they are for instance hints of horizontal gene transfer(s) with bacteria from other phyla.


----------------------------------------
Calculation of pairwise protein identity
----------------------------------------

Pairwise protein sequence identities reported on ChlamDB were calculated based on the multiple sequence alignments of orthologous groups (see `Orthogroup multiple sequence alignments & phylogenies`_). As the identity between two sequences can be calculated in different ways (see `this blog post discussing`_ ), we chose to calculate it by excluding gaps and calculating the identity based on aligned positions only:

* number of matches / ( number of matches +  number of mismatches)


------------------------------
Circular genome plots (Circos)
------------------------------

Circular genome plots are generated dynamically with Circos_. These plots are not generated based on the alignment of DNA sequences but based on orthology data (see previous paragraphs). The two outer gray circles show the location of `open reading frames`_ (ORFs) encoded on the leading and lagging strand of the reference genome (**Figure 2 A and B**). The red/blue inner circles show the conservation of each protein encoding ORF in one or multiple other genomes. Identity values were pre-calculated from protein alignments (see previous paragraphs). If the compared genomes encode more than one ortholog, the highest identity is used to color the region.

.. figure:: img/circos_method.png
    :figclass: align-center
    :width: 100%

    **Figure 2**: Example of circular genome plot. **A)** Compared genomes are ordered based on the median protein identity with the reference genome. **B)** Zoom showing the detail of a genomic region. All circular plots are interactive and users can click on any ORF to access the detailed annotation page of the corresponding protein.


------------------
Species phylogeny
------------------

The reference phylogeny is reconstructed with FastTree_ (default parameters, JTT+CAT model) based on the concatenated alignment of all single copy orthogroups conserved in all genomes.


-----------------------------------------------------------
Evaluation of the quality and completeness of draft genomes
-----------------------------------------------------------

The completeness and contamination of each genome is evaluated with checkM_ (v1.0.12). Completeness and contamination estimates are reported on the the main phylogeny page.


-------------------------------
Database and software versions
-------------------------------

Versions of softwares and databases used for a given zDB database are reported on the home page
of the webapplication.


.. _`biopython interface to Entrez`: https://biopython.org/DIST/docs/api/Bio.Entrez-module.html
.. _`boostrap support values`: https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1985.tb00420.x
.. _`Conserved Domain Database` : https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml
.. _`KEGG API`: https://www.kegg.jp/kegg/rest/keggapi.html
.. _`Kegg Orthologs`: https://www.genome.jp/kegg/ko.html
.. _`NCBI taxonomy ftp website`: ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
.. _`open reading frames`: https://en.wikipedia.org/wiki/Open_reading_frame
.. _`position-specific scoring matrix` : https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#CD_PSSM
.. _`this blog post discussing`: http://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
.. _`Virulence Factor Database`: http://www.mgc.ac.cn/VFs/
.. _AMRFinderPlus: https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/
.. _bioSQL: https://biosql.org/wiki/Main_Page
.. _BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi
.. _checkM: https://ecogenomics.github.io/CheckM/
.. _Circos: http://circos.ca/
.. _Diamond : https://www.nature.com/articles/nmeth.3176
.. _FastTree : http://www.microbesonline.org/fasttree/
.. _GenBank: https://www.ncbi.nlm.nih.gov/genbank/
.. _HMMER: http://hmmer.org/
.. _KOfam: https://www.genome.jp/tools/kofamkoala/
.. _KofamScan : ftp://ftp.genome.jp/pub/tools/kofamscan/
.. _MAFFT : https://mafft.cbrc.jp/alignment/software/
.. _MCL: https://micans.org/mcl/
.. _NextFlow: https://www.nextflow.io/
.. _OrthoFinder: https://github.com/davidemms/OrthoFinder
.. _Pfam: https://interpro-documentation.readthedocs.io/en/latest/pfam.html
.. _RefSeq: https://www.ncbi.nlm.nih.gov/refseq/
.. _rpsblast: https://www.ncbi.nlm.nih.gov/books/NBK279690/
.. _SwissProt: https://www.uniprot.org/
