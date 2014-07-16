# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up

- add orthomcl module 

module add SequenceAnalysis/OrthologyAnalysis/orthomclSoftware/2.0.9;

- add orthomcl path

echo 'export PATH=~/orthomcl/:$PATH' >> ~./bashrc

* Configuration
* Dependencies

Requires repository Utilitaires cloned and within PATH.

Requires blast module 

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

orthomcl.py -i file_1.faa file_2.faa ...

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact