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

v1.0.8 (february 2023):
- added support for dockers and conda
- several bugfixes

v1.0.5 (december 2022): 
- added a ```--resume``` option to ```zdb run```
- added a name field in the input csv file
- genbank files with other extensions are now accepted by the pipeline

## Installation

### zDB Installation

zDB can be installed from conda with the following command
```
conda install zDB -c metagenlab -c bioconda
```
For now, the project is hosted on our own conda channel. A bioconda package is also available, but is currently not up to date.

The analysis and the webserver can be run in either conda environments or in containers (both singularity and dockers are supported).

Running zDB in conda is likely the easiest way as it does not require to install either singularity or dockers, but comes with several drawbacks:
- django will be run in native mode, without nginx and gunicorn and should not be used to set up a web-facing database (it is fine for a local access)
- containers allow us to have a precise control of the environment where the webapp is run; it is less the case for conda environment. Despite our best care, running the webapp in conda might not work due to local differences.

Of note, zDB has been tested on singularity v3.8.3 and v3.8.4 but should work on more recent versions. 
If you opt to use singularity, it can be installed with the following command:
```
conda install singularity=3.8.4 -c conda-forge
```

### Install zDB from sources

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
```

### Too long, will not read

Here are a few examples of workflows with a test dataset available from the git repository.
To get the test dataset:
```
wget https://github.com/metagenlab/zDB/raw/master/test_dataset.tar.gz
tar xvf test_dataset.tar.gz
```

For a minimal database (assuming that singularity is installed):
```
conda install zdb -c metagenlab
zdb run --input=input.csv --name=simple_run # runs the analysis
zdb webapp --name=simple_run # Launches the webapp on simple run
```
The minimal database should take around 5 minutes to complete in a recent Desktop machine.

To do the same in conda environments:
```
conda install zdb -c metagenlab
zdb run --input=input.csv --name=simple_run_conda --conda # runs the analysis
zdb webapp --conda --name=simple_run_conda # Launches the webapp on the latest run
```

To have a more complete set of analyses (includes cog and pfam annotation):
```
conda install zdb -c metagenlab
zdb setup --pfam --cog --conda
zdb run --input=input.csv --name=more_complete_run --conda --cog --pfam # runs the analysis
zdb webapp --conda --name=more_complete_run # Launches the webapp on the latest run
```

## Setting up the reference databases

Depending on which analysis are to be run, reference databases will need to be downloaded and set up.
This is done using the ```zdb setup``` command.
You'll need to specifiy which databases are to be downloaded.
Of note, in minimal mode, zdb does not require any database to run.


The following databases can be downloaded:
```
--cog: downloads the CDD profiles used for COG annotations
--ko: downloads and setups the hmm profiles of the ko database
--pfam: downloads and setups up the hmm profiles of the PFAM protein domains
--swissprot: downloads and indexes the swissprot database
```

The database setup can be run either in singularity containers (by default), in conda environment (if the ```--conda``` flag is set) or in docker (if the ```--docker``` flag is set).

In addition, you can specify the directory where you want the databases to be installed with the ```--dir``` option: ```zdb setup --swissprot --dir=foobardir```.

Downloading the HMM files from the KEGG server can take a bit of time (but you'll only need to do this once).

Example commmand, setting up all the databases in the current directory, using conda environments:
```
zdb setup --pfam --swissprot --cog --ko --conda
```

## Running the analysis

Once you have the reference databases set up, the genomes ready, just run the ```zdb run``` command. The run command expects a csv file as input. The csv should look like:
```
name, file
,foo/bar.gbk
,baz/bazz.gbk
foobar,foobar/baz.gbff
```
The ```name``` column is optional and can be omitted from the input csv file. **By default, zdb will use the organism's name as defined in the genbank file to identify genomes in the web application**. Specifying a name for a genome will tell zdb to use that name instead of the organism name from the genbank file. This is practical when working with assembled genomes that haven't been named yet or when working with genomes of different strains of a same species. If the same name is used in different files, zdb will just add a numbering suffix to make the names unique.

Before launching the analysis, zdb will also check for the uniqueness of locus tags and generate new ones if necessary. This is usually not necessary for genomes downloaded from RefSeq or other databases, but if genomes were annotated with automated tools, name collisions might happen.

Several options are available and allow you to customize the run.

By default, the analysis are run in singularity containers, but you can change this by using the ```--conda``` or ```--docker``` flags to have them run in conda environments or docker containers, respectively. If singularity is enabled, the containers will have to be downloaded. By default, they are stored in the singularity folder of the current directory, but this can be changed using the ```--singularity_dir``` option. This might be useful if you want to share containers between analyses.

If the databases were set up, additional analyses can also be enabled with the ```--ko```, ```--cog```, ```--pfam``` and ```--swissprot``` flags. The directory (by default zdb_ref in the current directory) where the database were installed can be specified with the ```--ref_dir``` option.

Other options include:
```
--resume: wrapper for nextflow resume, allows to restart a run that crashed without redoing all the computations
--out: directory where the files necessary for the webapp will be stored
--input: CSV file containing the path to the genbank files to include in the analysis
--name: custom run name (defaults to the name given by nextflow). The latest completed run is also named latest.

--cpu: number of parallel processes allowed (default 8)
--mem: max memory usage allowed (default 8GB)
--singularity_dir: the directory where the singularity images are downloaded (default singularity in current directory)
```

The ```--name``` option is optional and can be used to replace nextflow's randomly generated run names by more meaningful ones.

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
```

The webserver can also be run in a conda environment by setting the ```--conda``` flag or in docker by setting the ```docker``` flag.
By default, it will be run in a singularity container.

Simply put, if the port 8080 is not in use, you can simply run the ```zdb webapp``` script without any parameters. zDB will launch the webapp on the last run of analysis.

Once the webserver has started, you'll be able to access the webpages with a web browser.

If you do not remember which runs are available, you can list them with the ```zdb list_runs``` command.

### Known issue

When trying to run several instances of the webapp on the same machine, you'll need to modify two files to avoid interferences between the instances.

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
