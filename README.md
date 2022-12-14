# zDB: comparative bacterial genomics made easy

zDB is designed to perform comparative genomics analysis and to integrate the results in a Django web-application.

Several analysis are currently supported, with more to come:

- COG annotation
- KEGG orthologs annotation
- PFAM domains annotation
- Swissprot homologs search
- (RefSeq homologs search): implemented, but significantly slows down the analysis. You'll also have to download and prepare the database for diamond search, as this was not included in the database setup script.

In addition, zDB performs orthology and phylogeny inference.
All the results are stored either in a SQLite database or directly as files and displayed in the web application.


## Changelog

v1.0.5 (december 2022): 
- added a ```--resume``` option to ```zdb run``
- added a name field in the input csv file
- genbank files with other extensions are now accepted by the pipeline

## Installation

zDB relies on singularity to run the analysis and the web server. Unfortunately, the singularity versions available in the bioconda channel are currently outdated and you'll need to install more recent ones from conda-forge. Run the following command:

```
conda install singularity=3.8.4 -c conda-forge
```
As of now, zDB has been tested with this version of singularity (and 3.8.3), but it should work on more recent releases.

Once this is done, zDB can be installed from conda with the following command
```
conda install zDB -c metagenlab -c bioconda
```

For now, the project is hosted on our own conda channel. A bioconda package is also available, but is currently not up to date.

### Install from sources
You can also install zdb directly from the github repository. This is particularly useful if you want to make modifications or if you want to have a direct access to Nextflow config file for a better control of the execution.

Check out the project or download and unpack a release, then edit this line of the bin/zdb bash script:
```
NEXTFLOW_DIR="${CONDA}/share/zdb-${VERSION}/"
```
and replace it by the directory where you downloaded the project (this should point to the directory where zdb's nextflow.config is located).
Add zdb's bin directory to PATH and voila, zdb should run smoothly.

## Overview

Several subcommands are available:
```
setup - download and prepare the reference databases
webapp - start the webapp
run - run the analysis pipeline
export - exports the results of a previous run in an archive
import - unpack an archive that was prepared with the export command in the current directory
       - so that the results can be used to start the webapp
list_runs - lists the completed runs available to start the website in a given directory
help - print this message
```

## Setting up the reference databases

Depending on which analysis are to be run, reference databases will need to be downloaded and set up. This is done using the ```zdb setup``` command.
You'll need to specifiy which databases are to be downloaded. The following options are available:
```
--cog: downloads the CDD profiles used for COG annotations
--ko: downloads and setups the hmm profiles of the ko database
--pfam: downloads and setups up the hmm profiles of the PFAM protein domains
--swissprot: downloads and indexes the swissprot database
```

In addition, you can specify the directory where you want the databases to be installed with the ```--dir``` option: ```zdb setup --swissprot --dir=foobardir```.

Downloading the HMM files from the KEGG server can take a bit of time (but you'll only need to do this once).

## Running the analysis

Easy. Once you have the reference databases set up, the genomes ready, just run the ```zdb run``` command.
Several options are available and allow you to customize the run:

```
--resume: wrapper for nextflow resume, allows to restart a run that crashed without redoing all the computations
--out: directory where the files necessary for the webapp will be stored
--input: CSV file containing the path to the genbank files to include in the analysis
--name: run name (defaults to the name given by nextflow). The latest completed run is also named latest.
--cog: perform cog annotation
--ko: perform ko (metabolism) annotation
--pfam: peform PFAM domain annotation
--swissprot: search for homologs in the swissprot database
--ref_dir: directory where the reference databases were setup up (default zdb_ref)
--cpu: number of parallel processes allowed (default 8)
--mem: max memory usage allowed (default 8GB)
--singularity_dir: the directory where the singularity images are downloaded (default singularity in current directory)
```
As the analysis are run in containers, nextflow will have to download the first time zDB is used. The containers will be stored in the singularity folder of the current directory. If you already downloaded the containers, you can point zDB to the location where they are located with the ````--singularity_dir`` to avoid new downloads.

The input CSV file should look like this:

```
name, file
,foo/bar.gbk
,baz/bazz.gbk
foobar,foobar/baz.gbff
```
and the command could look like this ```zdb run --input=my_genomes.csv --cpu=32 --pfam --out=my_outputdirectory```.

The ```name``` column is optional and can be omitted from the input csv file. By default, zdb will use the organism's name in the webapp. Specifying a name for a genome will tell zdb to use that name instead of the organism name from the genbank file. This is handy when working with assembled genomes that haven't been named yet.

Before launching the analysis, zdb will also check for the uniqueness of locus tags and generate new ones if necessary. This is usually not necessary for genomes downloaded from RefSeq or other databases, but if genomes were annotated with automated tools, name collisions might happen.

We noticed that undeterministic bugs sometimes happen when downloading a singularity container or when running long analysis. To resume the analysis when this happens, just add the ```--resume``` flag to the previous command. For example, if the run launched with the command  
```
zdb run --input=input.csv --ko --cog
```
crashed after 5 hours of analysis, you can resume the run to where it was before crashing with the command:
```
zdb run --input=input.csv --ko --cog --resume
```
This can also be used if you want to add another analysis to a previous run:
```
zdb run --input=input.csv --ko --cog --pfam --resume
```
in this case, only the pfam annotations will be performed as the other analysis have already completed.

## Starting the web server

Once the analysis is complete, the web application can be run with the ```zdb webapp``` command. The following options can be used:
```
--port=PORT_NUMBER      the port the application will be listening to, 8080 by default.
--name=RUN_NAME     when nextflow runs. If not specified, will default to the last successful run (latest).
--allowed_hosts=HOSTS   the name of the host or the ip address of the server. If not specified, will default to the ip addresses of the current host.
--singularity_dir: the directory where the singularity images are downloaded (default singularity in current directory)
```
Simply put, if the port 8080 is not in use, that you want to access the results of your latest run and will only access the web application locally, you can simply run the ```zdb webapp``` script without any parameters.

Once the webserver has started, you'll be able to access the webpages with a web-browser.

If you do not remember which runs are available, you can list them with the ```zdb list_runs``` command.

### Know issue

When trying to run several instances of the webapp on the same machine, you'll need to modify two files to avoid interferences between the instance.

In zdb/gunicorn/gunicorn.py:
```
bind = "0.0.0.0:8000"
```
modify the 8000 number to another, yet unused port.
Same in the zdb/nginx/nginx.config file:
```
proxy_pass http://localhost:8000;
```
Modify the 8000 to the same number you attributed to the port number of gunicorn.

## Importing and exporting results

As the analysis may be run on a server or on an HPC cluster, the results may need to be exported to start a web application on a different machine.

This can be done with the ```zdb export``` command with a run name as parameter. This will create a compressed archive containing all the necessary results. The archive can then be transferred to a different machine and unpacked, either manually or with the ```zdb import``` command.
The web server can then be started as if the analysis had been run locally.

## Bugs and feature requests
Suggestion and bug reports are very welcome [here](https://github.com/metagenlab/annotation_pipeline_nextflow/issues).

We already have several idea to improve the tool:
- make it possible to add or remove genomes in an existing database
- use d3.js to draw interactive phylogenetic trees
- add new annotations, in particular, we've already received some requests for antibiotic resistance and virulence

But we're definitely open for suggestions and contributions.
