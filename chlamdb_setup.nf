#!/usr/bin/env nextflow
/*
 * Author:
 * - Trestan Pillonel <trestan.pillonel@gmail.com>
 *
 */

log.info params.input
params.databases_dir = "$PWD/databases"
params.setup_COG = true
params.setup_enzyme = true

log.info "====================================="
log.info "input                  : ${params.input}"



Channel.from([["cognames2003-2014.tab", "ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab"],
              ["cog2003-2014.csv","ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv"],
              ["fun2003-2014.tab", "ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab"]])
              .set { cog_urls }

Channel.fromPath("${params.databases_dir}/annotation/COG/blast_COG.tab")
              .into { cog_results }


process download_COG {

  publishDir 'chlamdb_setup/COG_tables', mode: 'copy', overwrite: true
  echo true

  when:
  params.setup_COG == true

  input:
  set val (table_name), val (table_url)  from cog_urls

  output:
  file("${table_name}") into cog_tables

  script:
  """
  echo ${table_url}
  curl -L ${table_url} > ${table_name}
  """
}

process mysql_setup_COG_tables {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.setup_COG == true

  input:
  file cog_tables from cog_tables.collect()

  output:
  file("mysql_COG_setup.log") into mysql_COG_setup

  script:
  """
  echo ${cog_tables}
  chlamdb-setup-COG.py -i cognames2003-2014.tab -c cog2003-2014.csv -f fun2003-2014.tab > mysql_COG_setup.log
  """
}


process mysql_setup_enzyme_KEGG_tables {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.setup_enzyme == true

  output:
  file("mysql_enzyme_setup.log") into mysql_enzyme_setup

  script:
  """
  chlamdb-setup-enzyme-kegg.py -u > mysql_enzyme_setup.log
  """
}

process mysql_setup_biosql_db {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.setup_biosql == true

  output:
  file("mysql_biosql_setup.log") into mysql_biosql_setup

  script:
  """
  #!/usr/bin/env python

  import urllib.request
  import sys
  import MySQLdb
  import os

  log = open("mysql_biosql_setup.log", "w")

  url_mysql_biosql_scheme = 'https://raw.githubusercontent.com/biosql/biosql/master/sql/biosqldb-mysql.sql'

  log.write('Downloading Biosql scheme from %s ...\n' % url_mysql_biosql_scheme)
  request = urllib.request.Request(url_mysql_biosql_scheme)
  page = urllib.request.urlopen(request)

  log.write("Creating mysql database...\n")

  sqlpsw = os.environ['SQLPSW']
  conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                         user="root", # your username
                         passwd=sqlpsw) # name of the data base
  cursor = conn.cursor()
  sql_db = 'CREATE DATABASE IF NOT EXISTS biosqldb;
  cursor.execute(sql_db,)
  conn.commit()
  cursor.execute("use biosqldb;",)

  log.write("Importing Biosql schema...\n")
  conn.executescript(page.read().decode('unicode-escape'))

  conn.commit()
  log.close()
  """
}


workflow.onComplete {
  // Display complete message
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  mail = [ to: 'trestan.pillonel@gmail.com',
           subject: 'Annotation Pipeline - DONE',
           body: 'SUCCESS!' ]


}

workflow.onError {
  // Display error message
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
