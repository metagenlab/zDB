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

1. The first section 

2. The second section defines where the reference databases will be stored on disk. This is used both during the analysis and when setting up the reference databases. Unless you already have some of the reference databases installed, you won't have to modify this section.

3. The third section is where you define which analysis to run and the parameters for the tools that need them. By default, all analysis are enabled (set to true), except the search for refseq homologs. All analysis can be set to false in case the results are not relevant or if you don't want to install the reference database.

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

## Running the analysis

## Starting the web server

A default config file for nextflow is included in the repository and already sets up the 
