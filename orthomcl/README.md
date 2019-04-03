# README #

This README would normally document whatever steps are necessary to get your application up and running.
###TODO###

- change formatdb to new equivalent

### What is this repository for? ###

* Quick summary

Scripts to perform a clustering of amino acid sequences into groups of orthologous proteins.

* Version

V1.0

* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up

add orthomcl module 

```
#!bash

module add SequenceAnalysis/OrthologyAnalysis/orthomclSoftware/2.0.9;
```


add orthomcl path

```
#!bash

echo 'export PATH=~/bin/orthomcl/:$PATH' >> ~/.bashrc
```


* Configuration
* Dependencies

Requires repository utilitaires cloned and within PATH.
Add utilitaires to path

```
#!bash

echo 'export PATH=~/bin/utilitaires/:$PATH' >> ~/.bashrc
```

Requires blast module 

```
#!bash
module add Blast/blast/2.2.26;
module add Blast/ncbi-blast/2.2.28+;
```


1. scripts python:

- generate_bsub_file 

- splitFasta

- shell_command

## Database configuration ##


```
#!python
dbVendor=mysql
dbConnectString=dbi:mysql:orthomcl;host=sib-sql04.vital-it.ch
dbLogin=orthomcl
dbPassword=orthomcl
similarSequencesTable=similar_table
orthologTable=ortho_table_test
orthologTable=ortholog
inParalogTable=in_paralog
coOrthologTable=co_ortholog
interTaxonMatchView=inter_taxon_match
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE
```

* How to run tests
* Deployment instructions


## BASIC USAGE ##


```
#!python

orthomcl.py -i file_1.faa file_2.faa ...
```


### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin

claire.bertelli@chuv.ch
trestan.pillonel@chuv.ch

* Other community or team contact