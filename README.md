# zDB: comparative bacterial genomics made easy

zDB is a tool designed to allow 

Several annotations can be performed, depending on the settings:



## Installation

[Nextflow](https://www.nextflow.io/) and [singularity](https://sylabs.io/guides/3.5/user-guide/introduction.html) will need to be installed before running zDB (and of course, git).
Once both tools have been installed, checkout zDB with the following command:

```
git clone git@github.com:metagenlab/annotation_pipeline_nextflow.git
```

## Config file

The ```nextflow.config``` file is separated in several sections:

1. The input section. The annotation pipeline expects a simple tsv file as input (actually, could also be a csv file, as separators are currently not necessary). The file should list the genbank files to be included in the analysis under the ```file``` header. For example:

```
file
foo.gbk
bar.gbk
pof/baz.gbk
```

2. The second section defines where the reference databases will be stored on disk. This is used both during the analysis and when setting up the reference databases. Unless you already have some of the reference databases installed, you won't have to modify this section.

3. The third section is where you define which analysis to run and the parameters for the tools that need them. By default, all analysis are enabled (set to true), except the search for refseq homologs. All analysis can be set to false in case the results are not relevant or if you don't want to install the reference database. The core_missing parameter may be useful to build the species phylogeny in datasets that include incomplete genomes. If the parameter is set to 0, only the core genes present in all genomes will be taken into account to build the phylogeny.

4. The internals lists all the parameters necessary for pipeline to run. Modify those at your own risk and perils.

5. The last section allows you to specify the resources allocated to the analysis. You can limit CPU or memory usage by setting different values for the cpus or memory options.


## Setting up the reference databases

Depending on which analysis are to be run, reference databases will need to be downloaded and set up. 
This is done by running the db_setup.nf script with nextflow:

```
nextflow run db_setup.nf
```

The script will download the reference databases needed for the analysis set to true in the nextflow config file.
For instance, if only the swissprot homologs and COG annotation are set to true in the config file, only the swissprot and COG reference databases will be downloaded and prepared.

Be warned that setting up the KEGG orthologs can take a bit of time (read a few hours), as zDB needs to download the definition of all KEGG orthologs and modules.

## Running the analysis

Easy. Once you have the reference databases set up, the genomes ready and are happy with the nextflow.config file, just run the 
```
nextflow run annotation_pipeline.nf
```
command. You can also name the run with the ```--name=``` parameter, which may be useful in case you plan on running several analysis.

## Starting the web server

A default config file for nextflow is included in the repository and already sets up the 
