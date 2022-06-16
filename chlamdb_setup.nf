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
params.setup_biosql = true
params.setup_linear_taxonomy = false


log.info "====================================="
log.info "input                  : ${params.input}"


//
// setup reference tables
//



process mysql_setup_biosql_db {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.setup_biosql == true

  output:
  file("mysql_biosql_setup.log") into sql_db

  script:
  """
  chlamdb-setup-sqldb.py -d ${params.db_name} > mysql_biosql_setup.log
  """
}


process create_biodb {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.setup_biosql == true

  input:
  file sql_db

  output:
  val "biosqldb" into biosqldb

  script:
  """
#!/usr/bin/env python
from chlamdb.biosqldb import manipulate_biosqldb
manipulate_biosqldb.create_new_biodb("${params.db_name}")
  """
}

biosqldb.into{sql_db_to_COG
            sql_db_to_KEGG
            sql_db_to_taxonomy
            sql_db_to_interpro
            sql_db_to_annotations
            sql_db_to_TCDB
            } 


process download_COG {

  publishDir 'chlamdb_setup/COG_tables', mode: 'copy', overwrite: true
  echo true

  when:
  params.setup_COG == true

  input:
  file sql_db_to_COG


  output:
  file("cognames2003-2014.tab") into cog_names
  file("cog2003-2014.csv") into cog_list
  file("fun2003-2014.tab") into cog_fun

  script:
  """
  curl -L ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab > cognames2003-2014.tab
  curl -L ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv > cog2003-2014.csv
  curl -L ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab > fun2003-2014.tab
  """
}

process mysql_setup_COG_tables {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.setup_COG == true

  input:
  file("cognames2003-2014.tab") from cog_names
  file("cog2003-2014.csv") from cog_list
  file("fun2003-2014.tab") from cog_fun

  output:
  file("mysql_COG_setup.log") into mysql_COG_setup

  script:
  """
  chlamdb-setup-COG.py -b ${params.db_name} -i cognames2003-2014.tab -f fun2003-2014.tab > mysql_COG_setup.log
  """
}


process mysql_setup_enzyme_KEGG_tables {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73 beautifulsoup4=4.7.1 lxml=4.3.3'

  when:
  params.setup_enzyme == true

  input:
  file sql_db_to_KEGG

  output:
  file("mysql_enzyme_setup.log2") into mysql_enzyme_setup

  script:
  """
  echo ok
  chlamdb-setup-enzyme-kegg.py -d ${params.db_name} > mysql_enzyme_setup.log2
  """
}



process mysql_setup_linear_taxonomy {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.setup_linear_taxonomy == true

  input:
  file sql_db_to_taxonomy

  output:
  file("setup_linear_taxonomy.log") into mysql_linear_taxonomy_setup

  script:
  """
  chlamdb-setup-linear-taxonomy.py -d ${params.db_name} -s ${params.db_directory}/ncbi-taxonomy/linear_taxonomy.db -l setup_linear_taxonomy.log
  """
}

process mysql_setup_interpro {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  input:
  file sql_db_to_interpro


  when:
  params.setup_interpro == true

  output:
  file("mysql_setup_interpro.log") into mysql_setup_interpro

  script:
  """
  chlamdb-setup-interpro.py -d ${params.db_name} -v ${params.interpro_version} -d ${params.db_name} > mysql_setup_interpro.log
  """
}

process mysql_setup_TCDB {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  input:
  file sql_db_to_interpro


  when:
  params.setup_tcdb == true

  output:
  file("setup_TCDB.log") into mysql_setup_TCDB

  script:
  """
  chlamdb-setup-TCDB.py -d ${params.db_name} -l setup_TCDB.log
  """
}

//
// load annotations
//

process mysql_load_genbank {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  sql_db_to_annotations

  output:
  file("mysql_load_genbank.log") into mysql_load_genbank

  script:
  """
  chlamdb-load-gbk.py -g ${params.execution_dir}/data/gbk_edited/*gbk -d ${params.db_name} >> mysql_load_genbank.log
  """
}


process mysql_load_orthofinder {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file mysql_load_genbank

  output:
  file("mysql_load_orthogroups.log") into mysql_load_orthogroups

  script:
  """
  chlamdb-load-orthofinder.py -m ${params.execution_dir}/orthology/Results_output/WorkingDirectory/OrthoFinder/Results_output/Orthogroups/Orthogroups.txt -d ${params.db_name} > mysql_load_orthogroups.log
  """
}

process mysql_load_alignments {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file mysql_load_orthogroups

  output:
  file("mysql_load_alignments.log") into mysql_load_alignments

  script:
  """
  # number of cpu
  chlamdb-load-alignments.py -a ${params.execution_dir}/orthology/orthogroups_alignments/*faa -d ${params.db_name} -c 20 > mysql_load_alignments.log
  """
}

process setup_orthology_matrix {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file mysql_load_alignments

  output:
  file("orthology_mat.log") into mysql_orthology_mat
  file("orthology_mat_plas.log") into mysql_orthology_mat_plas

  script:
  """
  chlamdb-setup-comparative-tables.py -d ${params.db_name} -o >> orthology_mat.log
  chlamdb-setup-comparative-tables.py -d ${params.db_name} -o -a >> orthology_mat_plas.log
  """
}

process consenus_orthogroup_annot {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file mysql_orthology_mat

  output:
  file("orthogroups_annot.log") into orthogroups_annot

  script:
  """
  chlamdb-get-consensus-orthogroup-annotation.py -d ${params.db_name} -g >> orthogroups_annot.log
  """
}

process consensus_orthogroup_annot {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file orthogroups_annot

  output:
  file("old_locus_corresp.log") into old_locus_corresp

  script:
  """
  chlamdb-setup-old_locus-table.py -d ${params.db_name} >> old_locus_corresp.log
  """
}

process reference_phylo {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file orthogroups_annot

  output:
  file("reference_phylo.log") into reference_phylo

  script:
  """
  chlamdb-load-reference-phylogeny.py -r ${params.execution_dir}/orthology/core_alignment_and_phylogeny/core_genome_phylogeny.nwk -d ${params.db_name} -g ${params.execution_dir}/data/gbk_edited/*gbk >> reference_phylo.log
  """
}

process genomes_statistics {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file reference_phylo

  output:
  file("genomes_statistics.log") into genomes_statistics

  script:
  """
  chlamdb-setup-genomes-statistics.py -d ${params.db_name} >> genomes_statistics.log
  """
}

process orthogroup_phylo {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file reference_phylo

  output:
  file("orthogroup_phylo.log") into orthogroup_phylo

  script:
  """
  chlamdb-load-phylogenies.py -d ${params.db_name} -t ${params.execution_dir}/orthology/orthogroups_phylogenies_fasttree/*nwk >> orthogroup_phylo.log
  """
}

process inteproscan {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file orthogroup_phylo

  output:
  file("inteproscan.log") into inteproscan

  script:
  """
  chlamdb-load-interproscan.py -d ${params.db_name} -u ${params.execution_dir}/data/nr_mapping.tab -i  ${params.execution_dir}/annotation/interproscan/*tsv >> inteproscan.log
  """
}

process inteproscan_legacy {

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file inteproscan

  output:
  file("inteproscan_leg.log") into inteproscan_leg

  script:
  """
  chlamdb-load-interproscan.py -d ${params.db_name} -u ${params.execution_dir}/data/nr_mapping.tab -i  ${params.execution_dir}/annotation/interproscan/*tsv -l >> inteproscan_leg.log
  """
}

process inteproscan_TM_SP {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file inteproscan_leg

  output:
  file("inteproscan_TM_SP.log") into inteproscan_TM_SP

  script:
  """
  chlamdb-load-interproscan.py -d ${params.db_name} -u ${params.execution_dir}/data/nr_mapping.tab -a >> inteproscan_TM_SP.log
  """
}

process locus2hash {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file inteproscan_TM_SP

  output:
  file("locus2hash.log") into locus2hash

  script:
  """
  chlamdb-load-hash2locus.py -d ${params.db_name} -u ${params.execution_dir}/data/nr_mapping.tab >> locus2hash.log
  """
}

locus2hash.into {
  hash_to_pfam
  hash_to_uniref
  hash_to_swissprot
  hash_to_BBH_phylo
  hash_to_paperblast
  hash_to_protparams
  hash_to_GC_table
  hash_to_synteny
  hash_to_phyloprofile
  hash_to_synonyms
  hash_to_psortb
  hash_to_pdb
  hash_to_string_pmid
  hash_to_T3_effectors
}



process comparative_tables_pfam {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_pfam == true

  input:
  file hash_to_pfam

  output:
  file("comparative_tables_pfam.log") into comparative_tables_pfam
  file("comparative_tables_pfam_plasmids.log") into comparative_tables_pfam_plasmids

  script:
  """
  chlamdb-setup-comparative-tables.py -d ${params.db_name} -p >> comparative_tables_pfam.log
  # distinguish chromosome from plasmids
  chlamdb-setup-comparative-tables.py -d ${params.db_name} -p -a >> comparative_tables_pfam_plasmids.log
  """
}

process comparative_tables_interpro {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file comparative_tables_pfam

  output:
  file("comparative_tables_interpro.log") into comparative_tables_interpro
  file("comparative_tables_interpro_plasmids.log") into comparative_tables_interpro_plasmids

  script:
  """
  chlamdb-setup-comparative-tables.py -d ${params.db_name} -i >> comparative_tables_interpro.log
  # distinguish chromosome from plasmids
  chlamdb-setup-comparative-tables.py -d ${params.db_name} -i -a >> comparative_tables_interpro_plasmids.log
  """
}

process consensus_annot_og {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_annotations == true

  input:
  file locus2hash

  output:
  file("consensus_annot_og.log") into consensus_annot_og

  script:
  """
  chlamdb-get-consensus-orthogroup-annotation.py -d ${params.db_name} -i  >> consensus_annot_og.log
  """
}


process load_COG {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_COG == true

  input:
  file mysql_COG_setup

  output:
  file("load_COG.log") into load_COG

  script:
  """
  # chlamdb-load-COG.py [-h][-u HASH2LOCUS_TAG] [-cc COG2CDD] [-cl CDD2LENGTH]
  chlamdb-load-COG.py -d ${params.db_name} -i ${params.execution_dir}/annotation/COG/blast_COG.tab -u ${params.execution_dir}/data/nr_mapping.tab -cl ${params.cdd_length} -cc ${params.cdd_cog} >> load_COG.log
  """
}

process comparative_tables_COG {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_COG == true

  input:
  file load_COG

  output:
  file("comparative_tables_COG.log") into comparative_tables_COG
  file("comparative_tables_COG_plasmids.log") into comparative_tables_COG_plas

  script:
  """
  chlamdb-setup-comparative-tables.py -d ${params.db_name} -c  >> comparative_tables_COG.log
  chlamdb-setup-comparative-tables.py -d ${params.db_name} -c -a >> comparative_tables_COG_plasmids.log
  """
}

process consensus_annot_COG {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_COG == true

  input:
  file comparative_tables_COG

  output:
  file("consensus_annot_COG.log") into consensus_annot_COG

  script:
  """
  chlamdb-get-consensus-orthogroup-annotation.py -d ${params.db_name} -c >> consensus_annot_COG.log
  """
}

process load_KEGG {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_KEGG == true

  input:
  file mysql_enzyme_setup

  output:
  file("load_KEGG.log") into load_KEGG

  script:
  """
  chlamdb-load-KO.py -k ${params.execution_dir}/annotation/KO/chunk*.tab -d ${params.db_name} -c ${params.execution_dir}/data/nr_mapping.tab > load_KEGG.log
  """
}

process comparative_tables_KEGG {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_KEGG == true

  input:
  file load_KEGG

  output:
  file("comparative_tables_KEGG.log") into comparative_tables_KEGG
  file("comparative_tables_KEGG_plasmids.log") into comparative_tables_KEGG_plas

  script:
  """
  chlamdb-setup-comparative-tables.py -d ${params.db_name} -k  >> comparative_tables_KEGG.log
  chlamdb-setup-comparative-tables.py -d ${params.db_name} -k -a >> comparative_tables_KEGG_plasmids.log
  """
}

process consensus_annot_KEGG {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_KEGG == true

  input:
  file comparative_tables_KEGG

  output:
  file("consensus_annot_KEGG.log") into consensus_annot_KEGG

  script:
  """
  chlamdb-get-consensus-orthogroup-annotation.py -d ${params.db_name} -k >> consensus_annot_KEGG.log
  """
}

process load_TCDB {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_TCDB == true

  input:
  file mysql_setup_TCDB

  output:
  file("load_TCDB.log") into load_TCDB

  script:
  """
  chlamdb-load-TCDB.py -t GBLAST_HTML -b FASTA DATABASE  -k ${params.execution_dir}/annotation/KO/chunk*.tab -d ${params.db_name} -u ${params.execution_dir}/data/nr_mapping.tab > load_TCDB.log
  """
}


process load_uniref_hits {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_blast_uniref100 == true

  input:
  file hash_to_uniref

  output:
  file("load_uniref_hits.log") into load_uniref_hits

  script:
  """
  chlamdb-load-refseq-homology-search.py -t -f -d ${params.db_name} -u ${params.execution_dir}/data/nr_mapping.tab -lt ${params.linear_taxonomy} -rd  ${params.execution_dir}/annotation/diamond_uniref/diamond_uniref.db -ud ${params.uniref_db} > load_uniref_hits.log
  """
}

process load_swissprot_hits {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_blast_swissprot == true

  input:
  file hash_to_swissprot

  output:
  file("load_swissprot_hits.log") into load_swissprot_hits

  script:
  """
  chlamdb-load-swissprot-homology-search.py -l -p 8 -t -f -d ${params.db_name} -u ${params.execution_dir}/data/nr_mapping.tab -i ${params.execution_dir}/annotation/blast_swissprot/*tab > load_swissprot_hits.log
  """
}

process load_BBH_phylo {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_BBH_phylo == true

  input:
  file hash_to_BBH_phylo

  output:
  file("load_BBH_phylo.log") into load_BBH_phylo

  script:
  """
  chlamdb-load-phylogenies-BBH.py -d ${params.db_name} -t ${params.execution_dir}/orthology/orthogroups_uniref_diamond_BBH_phylogenies/*nwk > load_BBH_phylo.log
  """
}

process load_paperblast {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.load_paperblast == true

  input:
  file hash_to_paperblast

  output:
  file("load_paperblast.log") into load_paperblast

  script:
  """
  chlamdb-load-paperblast.py -d ${params.db_name} -p ${params.paperblast_sqlite} -t ${params.execution_dir}/annotation/paperblast_mapping/*tab -qf ${params.execution_dir}/data/nr.faa -df ${params.paperblast_fasta} -c ${params.execution_dir}/data/nr_mapping.tab > load_paperblast.log
  """
}

process protparams {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.protparams == true

  input:
  file hash_to_protparams

  output:
  file("protparams.log") into protparams

  script:
  """
  chlamdb-get-prot-params.py -d ${params.db_name} > protparams.log
  """
}


process GC_table {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.gc_table == true

  input:
  file hash_to_GC_table

  output:
  file("GC_table.log") into GC_table

  script:
  """
  chlamdb-setup-gc-content-tables.py -d ${params.db_name} > GC_table.log
  """
}

process synteny {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.synteny == true

  input:
  file hash_to_synteny

  output:
  file("synteny.log") into synteny

  script:
  """
  chlamdb-find-conserved-neighborhood.py -d ${params.db_name} > synteny.log
  """
}

process phyloprofile {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.phyloprofile == true

  input:
  file hash_to_phyloprofile

  output:
  file("phyloprofile.log") into phyloprofile

  script:
  """
  # both euclidian and jaccard dist
  chlamdb-get-orthogroup-profiles-euclidian-dist.py -d ${params.db_name} -c 8 -j -e  > phyloprofile.log
  """
}

process crossrefs {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.add_cross_references == true

  input:
  file hash_to_synonyms

  output:
  file("crossrefs.log") into crossrefs

  script:
  """
  chlamdb-load-cross-references.py -d ${params.db_name} -c ${params.execution_dir}/data/nr_mapping.tab -r ${params.execution_dir}/data/refseq_corresp/refseq_corresp.tab -i ${params.execution_dir}/annotation/uniparc_mapping/uniparc_crossreferences.tab > crossrefs.log
  """
}

process load_psortb {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.psortb == true

  input:
  file hash_to_psortb

  output:
  file("load_psortb.log") into load_psortb

  script:
  """
  chlamdb-load-psortdb.py -d ${params.db_name}  -c ${params.execution_dir}/data/nr_mapping.tab -t ${params.execution_dir}/annotation/psortb/*txt > load_psortb.log
  """
}


process load_pdb {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.pdb == true

  input:
  file hash_to_pdb

  output:
  file("load_pdb.log") into load_pdb

  script:
  """
  chlamdb-load-psortdb.py -d ${params.db_name}  -c ${params.execution_dir}/data/nr_mapping.tab -t ${params.execution_dir}/annotation/pdb_mapping/chunk*tab > load_pdb.log
  """
}

process load_T3_effectors {
  // # update TM et SP columns from legacy `ortho_detail` table

  publishDir 'chlamdb_setup/logs', mode: 'copy', overwrite: true
  echo true
  //conda 'mysqlclient=1.3.10 biopython=1.73'

  when:
  params.t3_effectors == true

  input:
  file hash_to_T3_effectors

  output:
  file("t3_effectors_BPBAac.log") into t3_effectors_BPBAac
  file("t3_effectors_DeepT3.log") into t3_effectors_DeepT3
  file("t3_effectors_effectiveT3.log") into t3_effectors_effectiveT3
  file("t3_effectors_T3_MM.log") into t3_effectors_T3_MM

  script:
  """
  chlamdb-load-effector-pred.py -d ${params.db_name}  -c ${params.execution_dir}/data/nr_mapping.tab -t ${params.execution_dir}/annotation/T3SS_effectors/BPBAac_results.tab -b > t3_effectors_BPBAac.log
  chlamdb-load-effector-pred.py -d ${params.db_name}  -c ${params.execution_dir}/data/nr_mapping.tab -t ${params.execution_dir}/annotation/T3SS_effectors/DeepT3_results.tab -dt > t3_effectors_DeepT3.log
  chlamdb-load-effector-pred.py -d ${params.db_name}  -c ${params.execution_dir}/data/nr_mapping.tab -t ${params.execution_dir}/annotation/T3SS_effectors/effectiveT3_results.tab -eff > t3_effectors_effectiveT3.log
  chlamdb-load-effector-pred.py -d ${params.db_name}  -c ${params.execution_dir}/data/nr_mapping.tab -t ${params.execution_dir}/annotation/T3SS_effectors/T3_MM_results.tab -t3 > t3_effectors_T3_MM.log
  """
}

//  STRING links to litterature





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
