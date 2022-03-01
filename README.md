# zDB: comparative bacterial genomics made easy

# Installation

Only 

# Config file

# Setting up the reference databases

Depending on which analysis are to be run, reference databases will need to be downloaded and set up. 
This is done by running the db_setup.nf script with nextflow:

```
nextflow db_setup.nf
```

The script will download the reference databases needed for the analysis set to true in the nextflow config file.
For instance, if only the swissprot homologs and COG annotation are set to true in the config file, only the swissprot and COG reference databases will be downloaded and prepared.


# Running the analysis

# Starting the web server


A default config file for nextflow is included in the repository and already sets up the 
