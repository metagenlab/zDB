================
Website tutorial
================

.. warning:: 
    zDB is still actively being developped and this tutorial might not always reflect the latest state of the application. If you have any question or suggestion feel free to open an issue on the `zDB github repository`_ or contact us (trestan.pillonel@chuv.ch or niklaus.johner@chuv.ch).

This page attempts to provide the user with help to understand and perform the analyses offered in the webinterface.

In this tutorial a database composed of a collection of *Enterobacteriaceae* genomes is presented with a focus on the fliLMNOPQR operon and specifically on FliL. This operon is FlhDC-dependent and encodes important regulatory factors and structural components of the membrane-spanning basal body of the flagellum.


--------------------------------
HOME page
--------------------------------
On the HOME page, the user has an overview of the database, notabaly a summary of the number of genomes and annotations in the database, as well as a box showing the status of the DB, notably the the versions of the reference databases used for the various annotations.

It also provides descriptions of the available analyses with direct links to the analysis views. The same links are also found in the left menu to facilitate the navigation.

Note that it is necessary to re-run the pipeline to incude the new analyses if wished (you can use the `--resume` argument - see documentation)


--------------------------------
Genome table and phylogeny
--------------------------------
The genome table is the first entry point evaluate the content of the database, getting details about the contigs and the loci identified on each genome (the clickable locus tags redirect to the `Protein annotation view`_ page).
The phylogeny, built on concatenated single copy orthologs with FastTree, shows the evolutionary relationships between the genomes given in input and essential data for the quality assessment of the given sequences.


--------------------------------
Homology search - Blast
--------------------------------
Perform a blast search of one or more sequences of interest against one or more genomes of the database.

Either an amino-acid or a nucleotide sequence can be given as input.

Set up:

* the type of homology search according to the input file:
     - blastp, tblastn with an aa sequence
     - blast_ffn, blast_fna, blastx with a nt sequence
* e-value
* maximum number of hits to display
* target genome or all genomes

.. note::
    If the search is performed against all genomes, max number of hits should be set to 'all' to avoid losing high identity matches.


In the reported example the protein sequence of the genes of fliLMNOPQR operon extracted from the *Enterobacter soli* genome are blasted (blastp) against all genomes of the database (:numref:`fig_blastp`).

.. _fig_blastp:

.. figure:: img/blastp_search.svg
    :figclass: align-center
    :width: 100%

    Blast interface for homology search. Blastp of the fliLMNOPQR operon genes (target: all, max number of hits: all). In the 'blast input' box amminoacid sequences of all genes introduced by a header.


This analysis allows to identify whether any of these genes are present in the genomes and evaluate the quality of the alignment of each hit (:numref:`fig_blastp_res2` Result 1):

| **A**. visual identification of hits for fliL gene
| **B**. info table about the hits (genome, contig/locus_tag, alignment scores and identity - Note that the locus tags are clickable and linked to the `Protein annotation view`_),
| **C**. Alignment of the query and the sequence of a hit selected in table B.

Additionally, the generated annotated phylogeny facilitate the interpretation of their distribution and conservation along all the genomes. As shown in :numref:`fig_blastp_res2` Result 2, four genomes carry all the investigated genes, fourteen genomes do not carry them, while the remaining ones have an incomplete set.

 
.. figure:: img/blastp_result1.svg
    :figclass: align-center
    :width: 100%

.. _fig_blastp_res2:

.. figure:: img/blastp_result2.svg
    :figclass: align-center
    :width: 100%

    **Blastp results** . *Result 1*: Details and *Result 2*: Phylogenetic distribution.


.. hint::
    - If you are interest in a specific gene expected to be present in one of the genomes included in the database, you can either retrieve the sequence in a public database, such as SwissProt, or use the search bar in the left-side menu of the web interface. Type the gene name, and identify which loci are annotated with that gene, clicking on one of them the user can directly retrieve both the nucleotide and the amino acid sequence of the gene - see `Protein annotation view`_ page below.
    - Compare the genomic regions around a protein of interest in selected genomes accessing the 'MENU/Genome alignments/Plot region' analysis - see the `Genome alignments`_ page below.


--------------------------------
Annotations
--------------------------------
The exact set of annotations available for analysis depends on the settings used during the generation of the database - see `Running the analysis`_ for an extensive explanation.

It allows the user to compare several aspects of selected genomes and perform comparative analyses for each annotation type: a) Orthogroups are always included in the database. The following annotations are optional and therefore not always available: b) KEGG orthologs, c) COG cluster, d) PFAM domains, e) Virulence factors, and f) Antimicrobial resistance genes.

Before proceeding here a brief description of the mentioned annotations and the link to the corresponding reference databases:
    * **Kegg**: Kegg annotations refer to the Kyoto Encyclopedia of Genes and Genomes (KEGG_). The genome annotation is composed of two aspects: a) KO assignemnt (KO is the identifier given to a functional ortholog defined from experimentally characterized genes and protein in specific organism), b) KEGG mapping where each KO is stored in a PATHWAY or MODULE identified based on molecular networks. This database provides a highly curated and reliable description of the metabolic pathway of the annotated genomes.
    * **COG**: COG annotations refer to the database of Cluster of Orthologous Genes (COGs_). In this database each COG is assigned to a functional category including metabolic, signal transduction, repair and other pathways. This database allows an easy comparison of organisms based on their preference for certain pathways.
    * **Pfam**: Pfam annotations refer to the Pfam_ database used to identify protein families and domains. Due to the nature of proteins as combinations of fixed structure, this database is based on the idea that the identification of domains wihin proteins can provide insights about their function.
    * **VFs**: Virulence factor annotations refer to the VFDB_, a curated database of virulence factors of bacterial pathogens.
    * **AMRs**: AMR annotations are obtained with the AMRFinderPlus_ software, which returns not only antimicrobial resistance genes, but also virulence and stress resistance genes. AMRFinderPlus_ uses highly curated gene and HMM databases from NCBI.

.. note::
    The following example in **Fig. 3 refers to the Orthogroups analyses**, however the same the same analyses are available for KEGG, COG, PFAM domains, VFs or AMR genes (check the help paragraph entitled '*Additional plots for Kegg Orthologs and Cluster of Orthologous Groups (COGs)*' to discover the extra analyses available with some of these annotations).


Overview of Orthogroups analyses
================================

Orthogroups are identified with Orthofinder_, an accurate platform that cluster *set of genes that are descended from a single gene in the last common ancestor of all the species being considered* as reported in its publication_.
In the following example, the orthogroup content is compared between *Enterobacter soli, Enterobacter ausbriae, Enterobacter ludvigii, and Klebsiella variicola* genomes. 


Extraction Form
---------------
Identify those orthogroups uniformly present in a set of genomes of interest and, optionally absent in others. Flexibility can be given to include orthogroups that, although present in some of the selected genomes, are not uniformly present in all and are missing in some ('Missing data' parameter).

The results box contains:

* A summary of the selected settings for the comparative analysis (:numref:`fig_extract` 1a): the orthgroup of 4 genomes are compared, no orthogroup will be excluded if present in other genomes, orthogroup that are present in 3 out of the 4 selected genomes are also reported.
* A list of identified orthogroups, description and distribution in the selected genomes (:numref:`fig_extract`1b): clicking on a Orthogroup entry redirects the user to the *Orthogroup annotation summary* page.
* A list of locus tags per each orthogroup and genome (:numref:`fig_extract` 1c): clicking on a Orthogroup entry redirects the user to the `Protein annotation view`_ page.

.. _fig_extract:

.. figure:: img/OverviewOrt_r1_r2.svg
    :figclass: align-center

    Orthogroups comparison overview of *Enterobacter soli, Enterobacter ausbriae, Enterobacter ludvigii, and Klebsiella variicola*. Analysis 2, 3, and 4 are reported in Fig. 4; analysis 5 is reported in Fig. 5.


Venn diagram
------------
Select a maximum of 6 genomes to visualize the distribution of their Orthogroups. This representation (:numref:`fig_venn_size_heat` 2) simplifies the identification of similarity/dissimilarity of Orthogroups between a few genomes.


Compare Orthogroup size
-----------------------
Visualize the number of entries of each Orthogroup in common between a selected set of genomes. This representation (:numref:`fig_venn_size_heat` 3) higlights which orthogroups are enriched or poorly represented in the genomes of interest.


Whole proteome heatmaps
-----------------------
Heatmap (:numref:`fig_venn_size_heat` 4) of presence/absence of the pool of Orthogroups present in the selected genomes. Discover which Orthogroups are widely shared by a subset of interest and which genome differentiate from the others. Going over the plot with the mouse it displays the orthogroup name, the organism of interest and the nummber of hits associated to that Orthogroup.

.. _fig_venn_size_heat:

.. figure:: img/Ort_venSize_heat.svg
    :figclass: align-center

    Orhogroup comparison analyses of *Enterobacter soli, Enterobacter ausbriae, Enterobacter ludvigii, and Klebsiella variicola*.


Pan/Core genome plot
--------------------
Graphical representation (:numref:`fig_core_pan_ort` 5) of the pan- and core- genome of a subset of genomes or of the uploaded dataset (**Fig. 5**).

This analysis generates three plots that display the content and conservation of Orthologous groups in selected genomes of interest.

    * *Total number of Orthologous groups* (:numref:`fig_core_pan_ort` 5A): this plot shows the number of all Orthologous groups present in a set of genomes. If the green curve reaches a plateau we can talk about 'closed pangenome' since no new Orthogroups are carried by additional genomes, on the contrary if the increment of the curve grows when looking at other genomes we can talk about 'open pangenome'.
    * *Number of shared Orthologous groups* (:numref:`fig_core_pan_ort` 5B): The red curve represents the core Orthogroups shared by the genomes and it tends to decrease as much as the compared genomes are different.
    * *Distribution of Orthologous groups* (:numref:`fig_core_pan_ort` 5C): the blue curve represents the number of Orthologous groups present in exactly n genomes displayed in the x-axis. This representation is useful to appreciate how many Orthologous groups are present in the totality of the genomes of interest, for example, or the diversity brought by single genomes. For example, if tot-1 is low it means that there are no specific genomes that bring a unique Orthologous groups.

.. _fig_core_pan_ort:

.. figure:: img/Core_pan_Ort_three.svg
    :figclass: align-center

    Accumulation/rarefaction plots.


Additional plots for Kegg Orthologs and Cluster of Orthologous Groups (COGs)
============================================================================

The comparative analyses of Kegg orthologs and COGs come with additional plots:

Categories Barchart
-------------------

Barchart plot (:numref:`fig_cog_barchart` 1) of the distribution of the entries annotated with a COG/KEGG category of selected genomes. It allows the evaluation of potential increment or decrement of entries known to be relevant for a certain function in some genomes of interest (:numref:`fig_cog_barchart` 1).

Focusing on the COG 'Cell motility' category, we see that *Klebsiella variicola* has fewer annotations of that category than *Enterobacter soli, Enterobacter ausbriae*, and *Enterobacter ludvigii*.

.. _fig_cog_barchart:

.. figure:: img/COGs_overview_bar_o.svg
    :figclass: align-center

    COGs comparison page. Barchart for each COG category representing the number of entries identified in each genome. The 'Cell motility' category is highligthed in green to stress the differences between the four selected genomes. Analyses 2 and 3 are reported in :numref:`fig_cog_heatmap`.


Categories Heatmaps
-------------------
Heatmaps of the COGs (not available for Keggs) along all the genomes expressed as fequency (:numref:`fig_cog_heatmap` 2*) or number (:numref:`fig_cog_heatmap` 3*) of identified entries.

Here the focus is again on the COG 'Cell motility' category where it emerges that *Klebsiella variicola* has 67 loci annotated in this category that represents 1.29% of total number of its loci, while *Enterobacter soli* has more than the double of its loci annotated in this category, 2.76% of them.

.. _fig_cog_heatmap:

.. figure:: img/COGs_heatmaps_o.svg
    :figclass: align-center

    Heatmaps of presence/absence of entries annotated with each COG category expressed as counts (2) or as frequencies (3). In the green box, the 'Cell motility' category, in purple, the two genomes of interest.


--------------------------------
Genome alignments
--------------------------------
This set of analyses allow the user to align the genomes and check the conservation of specific regions of interest.

Two plots can be generated:
    * circos
    * Plot region


Circos
======
Genomes alignment visualized in an interactive circular layout. This plot can trigger the identification of differentially distributed genomic regions in the genomes of interest, the presence of potential plasmid(s), or the products of other HGT events when looking at the GC composition, for example.
Following the help box, it is possible to recognize which regions encode for genes or tRNA and evaluate the conservation of the sequence checking the identity percentages.

In :numref:`fig_cog_heatmap`A, *Enterobacter ausbriae, Enterobacter ludvigii, and Klebsiella variicola* are mapped against 'Enterobacter soli'. The genomes appears similar in terms of gene content, however *Enterobacter soli* carries a plasmid which is absent in the other genomes.
When the user clicks on a gene of interest the `Protein annotation view`_ page will be displayed and provide the user with all the information about function, distribution and conservation of this protein.

.. note::
    the regions present in one of the compared genomes but in the reference, will not be visualized. A new plot inverting the genome given as reference will give this info.


Plot region
===========
'Plot region' analysis allows the user to discover a specific genomic region of interest. It plots the genomic features located in the neighborhood of a provided target locus, it displays the conservation of the protein of interest and the genes present in the flanking region among selected genomes (max 20000 bp).

In :numref:`fig_cog_heatmap`B, the focus is on the fliL gene of the fliLMNOPQR operon in *Enterobacter soli, Enterobacter ausbriae, Enterobacter ludvigii, and Klebsiella variicola*.
The operon is highly conserved in the Enterobacter genomes, but absent in *Klebsiella variicola*, which is indeed not reported in the plot (:numref:`fig_cog_heatmap`B). (Note that the phylogeny obtained in *Homology search - Blast*, already highlight the lack of these genes in *Klebsiella variicola* ).

.. _fig_plot_region:

.. figure:: img/Plot_region_ENTAS_RS13815_fliL_Soli_o_vertical.svg
    :figclass: align-center
    :width: 100%

    **Figure 8.** A. Circos plot of four genomes of interest and B. focus on the genomics region (20000 bp) around fliL gene (fliLMNOPQR operon). The operon is conserved among Enterobacter soli, Enterobacter asburiae and Enterobacter ludwigii. In red the gene encoded in the locus tag provided, in green CDs, in black the pseudogenes, and in yellow rRNAs and tRNAs.


--------------------------------
Metabolism
--------------------------------
This section provides the user with a set of analyses useful to discover the metabolism of given genomes based on the KEGG Orthology database.
It relies on the functional orthologs of the KO database which are categorized in molecular interaction, reaction and relation networks, named *KEGG pathway maps*, and functional units of gene sets, named *Kegg modules* associated with metabolism.


Kegg maps
=========
With this analysis the **Kegg pathways** of a genome of interest can be discovered, which Kegg orthologs of the pathway are present and compare their distribution in the other genomes.
In the following example (:numref:`fig_metabo_kegg_maps`), the Kegg pathways present in the *Enterobacter Soli* genome are listed (235 pathways in total) and a heatmap of the Ko of the flagellar pathways is shown. In this page a direct link to the official Kegg page is provided to evaluate the state of composition of this Kegg map (in red the KOs present in *Enterobacter soli*.

.. _fig_metabo_kegg_maps:

.. figure:: img/Metab_kegg_maps_o.svg
    :figclass: align-center

    Metabolism/kegg maps analysis. Steps to identify the completeness of a Kegg pathway for a genome of interest. The flagellar assembly pathways of *Enterobacter soli* is shown.


Kegg modules
============
Discover the KO of Kegg modules, organized in categories and sub categories, of a genome of interest or a subset of them (:numref:`fig_metabo_kegg_modules`).
Three types of search are available:

| **Category heatmap**: discover a Kegg category of interest, such as Energy metabolism and get an overview of the presence/absence of the kegg modules part of this category in the whole set of genomes. KO entry M00175 refers to 'Nitrogen fixation, nitrogen --> ammonia and it is present only in a few genomes, and one of them is *Klebsiella variicola* (:numref:`fig_metabo_kegg_modules` A).
| **Sub category heatmap**: similar output than the 'Category heatmap' search, but considering subcategories - for example ATP synthesis.
| **Compare strains**: this search let the user focus on a selected set of genomes to compare all the Kegg modules carried by them and better appreciated their distribution within the genomes. In :numref:`fig_metabo_kegg_modules` B, the four genomes are compared.

.. _fig_metabo_kegg_modules:

.. figure:: img/Metab_kegg_modules_Energy_met_o.svg
    :figclass: align-center

    Metabolism/kegg module analysis. A 'Category heatmap' output, B: 'Compare strains' output.

.. note::
    Search 1 and 3 come with a link to the `Kegg module overview`_ page.


Kegg module overview page
=========================
This page is accessible clicking on the Kegg module entry from the 'Metabolism/Kegg module' analysis or from the 'Locus tag overview page'. It gives access to the list of Ko entries that form the Kegg module of interest, and provides an indication of the completeness of the Kegg module within the genomes of the database.

The reported example is based on the KO entries of the kegg module number M00049 which describes the Adenine ribonucleotide biosynthesis ( IMP => ADP,ATP), and it is part of the *Nucleotide metabolism* category and *Purine metabolism* subcategory. Four genes are required to have a complete module, and one of them can be one among a set of four redundant genes. Among the genomes of the dataset, all except three have a complete module.

.. _fig_kegg_overview_page:

.. figure:: img/kegg_overview_page_IMP_o.svg
    :figclass: align-center

    Phylogeny annotatedd with presence/absence of KO entries of kegg module M00049.


------------------------
Protein annotation view
------------------------
This page provides a complete overview of a selected locus of interest.
The annotations are automatically retrieved from the .gbk files given as input, while further annotations can be assigned with COG, KEGG, Pfam, Swissprot, and Refseq databases only upon request (Note that RefSeq annotations are highly computational- and time-demanding)

In the example reported (:numref:`fig_locus_tag_overview`), the page displays the locus tag ENTAS_RS13815 of *Enterobacter soli* annotated with the fliL gene. The following info can be retrieved from the 'Overview' page:

- 1. A summary of the locus tag name, its size, the gene name if annotated and gene product are reported.
- 2. The Orthologous group to which the locus tag is assigned, the number of homologs of that orthogroup, the number of genomes in which the orthogroup is present.
- 3. the genomic region around the locus tag of interest. This plot provides an interactive way to discover of the flanking region of the target.
- 4. Box with useful functional and metabolic annotations (adatpted to the requested annotations in the config file)

.. _fig_locus_tag_overview

.. figure:: img/Locus_tag_filL_overview_m_o.svg
    :figclass: align-center

    Locus tag overview page. Overview of the locus tag ENTAS_RS13815 of *Enterobacter soli* encoding fliL gene.

From the 'Overview' page further plots are accessible (:numref:`fig_locus_tag_plots`):
the phylogenetic distribution of the orthogroup of the locus tag (**A**),the homologs of which are reported in a phylogeny with a dedicated attention to the Pfam domains composing them (**D**). Additionally, SwissProt and RefSeq annotations are listed to further evaluate the best homologs according to their databases (**B** and **C**) and the best RefSeq hits are included in the homologs phylogeny (**E**).
These analyses better characterize the locus whether the other annotations are not consistent for example, to infer horizontal gene transfer occurences, and also to observe potential dissimilarities/similarities in terms of Pfam domains between members of the same orthogroup. 

.. _fig_locus_tag_plots

.. figure:: img/Locus_tag_filL_plots_m_o.svg
    :figclass: align-center

    **Figure 13: Locus tag page plots**. A: phylogenetic distribution of the orthogroup; B: Homologs of ENTAS_RS13815 locus tag identified in RefSeq; C: Homologs of ENTAS_RS13815 locus tag identified in SwissProt; D: Orthogroup phylogeny of group_2742 with Pfam domains annotation; E: Phylogeny of the orthogroup identified in the set of genomes plus the addition of the three best RefSeq hits of locus tag ENTAS_RS13815.


.. note::
    In the boxes with Kegg, COGs, and Pfam annotations, you will be redirected to their explanatory overview pages (3 ouputs, all similar, with link to external sources, occurences in proteins in the orthologous groups, then list of locus tags with that annotation in all the genomes of the database, phylogeny of the dataset annotated with the number of hits for that annotation and their distribution in the orthologous groups --- MAYBE PUT AN EXAMPLE OF THAT PAGE FOR ONE ANNOTATION  )


-----------------------------
Orthogroup annotation summary
-----------------------------
This page represents several overlaps with the `Protein annotation view`_ page, however this is focused on the orthogroup rather than on a single member and its homologs. Indeed, it may occur that the homologs of a locus tag are split within more orthogroups.
Of interest, in this page the alignment between the members of the orthogroup is available and amino acid substitutions can be easily observed (:numref:`fig_og_overview` A)

.. _fig_og_overview

.. figure:: img/Orthogroup_page_overview_align_m_o.svg
    :figclass: align-center

    Overview of orthogroup 2742 of fliL gene of *Enterobacter soli* and protein alignment of its members.


------------------------------
KO/COG/Pfam annotation summary
------------------------------
A summary page of each COG, Pfam, and Kegg entry is accessible in the web interface through the analysis in the ``Comparison`` section pages, through the `Protein annotation view`_ page and even from the ``Metabolism`` section pages.
Each page provides a complete overview of the investigated annotation within the database and it comes also with external links.

It is organized in three sections that can be visualized in :numref:`fig_pfam_overview` where Pfam domain PF03748 is reported:
    * **General**: It provides how many loci are characterized with that annotation combining the info with the Orthogroups classification.
    * **Protein list**: list of all locus tags with that annotation
    * **Profile**: phylogeny annotated with an heatmap of the entries with that annotation and their distribution into Orthogroups

.. _fig_pfam_overview

.. figure:: img/Pfam_overview_page_o.svg
    :figclass: align-center

    Overview of Pfam domain PF03748.


--------------------------------
Search bar
--------------------------------
The search bar at the top of the left-side menu recognizes the following entries:

=============================   =================
Name 	                        Example
=============================   =================
KO entry             	        K02415
COG entry                    	COG1580
COG name                        Glutamate-1-semialdehyde aminotransferase
Gene name 	                    fliL
AMR genes                       ampC
Vilurence factors               VFG049129
Gene product 	                flagellar basal body-associated protein FliL
Locus tag accession name 	    ENTAS_RS13815
Organism	                    Enterobacter soli
=============================   =================

It is built with Whoosh_ and it can take in input also combination of terms separated by AND/OR, for a more complex search, for example. 

.. _publication : https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0721-2
.. _Orthofinder : https://github.com/davidemms/OrthoFinder
.. _KEGG : https://www.genome.jp/kegg/ko.html
.. _COGs : https://www.ncbi.nlm.nih.gov/research/cog
.. _Pfam : http://pfam.xfam.org/
.. _Whoosh : https://whoosh.readthedocs.io/en/latest/index.html
.. _`zDB github repository`: https://github.com/metagenlab/zDB
.. _VFDB: http://www.mgc.ac.cn/VFs/
.. _AMRFinderPlus: https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/
.. _`Running the analysis`: https://zdb.readthedocs.io/en/nj-docs/include_readme.html#running-the-analysis
