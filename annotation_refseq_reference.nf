#!/usr/bin/env nextflow

log.info params.input
params.databases_dir = "$PWD/databases"

log.info "====================================="
log.info "Database folder        : ${params.databases_dir}"


def get_interpro_version(){
    def url = 'https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload';
    def page = url.toURL().text;
    def version = (page =~ /interproscan-(\d.\d+-\d+.\d)/);
    return version[0][0];
}


process download_refseq_table {

  publishDir 'databases/refseq/', mode: 'copy', overwrite: false

  output:
    file 'assembly_summary_refseq.txt' into refseq_table

  script:
  """
  curl -S ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt | tail -n +2 | sed 's/# assembly_accession/assembly_accession/' > assembly_summary_refseq.txt
  """
}


/* .flatten().map { it }.filter { (it.refseq_category ==~ "reference genome")} */
refseq_table.splitCsv(sep: '\t', header: true).filter { (it.refseq_category ==~ "reference genome|representative genome")}.set { representative_and_reference_genomes }

process download_reference_genomes {

  publishDir 'databases/refseq/genomes', mode: 'copy', overwrite: false

  input:
    each assembly from representative_and_reference_genomes

  output:
    file "${assembly.assembly_accession}.faa.gz" into assembly_faa_gz

  script:
  def asm_name = "${assembly.asm_name}".replace(" ", "_")
  def asm_last = "${assembly.ftp_path}".split('/').last()
  """
  echo LAST ${asm_last}
  echo ${assembly.ftp_path}/${asm_last}_protein.faa.gz -o ${assembly.assembly_accession}.faa.gz
  curl  ${assembly.ftp_path}/${asm_last}_protein.faa.gz -o ${assembly.assembly_accession}.faa.gz
  """
}

process extract_gz {

  publishDir 'databases/refseq/genomes', mode: 'copy', overwrite: false

  input:
    file assembly_faa_gz

  output:
    file "${assembly_faa_gz.baseName}" into assembly_faa

  script:
  """
  zcat ${assembly_faa_gz} > ${assembly_faa_gz.baseName}
  """
}


assembly_faa.flatten().map { it }.filter { (it.text =~ /(>)/).size() < 1000 }.set { small_genomes }

process execute_interproscan {

  publishDir 'annotation/interproscan_results', mode: 'copy', overwrite: false

  cpus 4

  input:
  each file(seq) from small_genomes

  output:
  file '*tsv' into interpro_tsv

  script:
  n = seq.name
  """
  interproscan.sh -cpu 4 -appl ProDom,SFLD,HAMAP,TIGRFAM,SUPERFAMILY,PRINTS,PIRSF,ProSiteProfiles,PfamA,SMART,CDD -f TSV -i ${n} -d . -T . -iprlookup
  """
}
