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

1. The first section defines where the reference databases will be stored on disk. This is used both during the analysis and when setting up the reference databases. Unless you already have some of the reference databases installed, you won't have to modify this section.
2. 

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
