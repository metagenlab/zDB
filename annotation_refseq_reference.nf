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
    file "${assembly.assembly_accession}.faa.gz" into assembly_faa

  script:
  def asm_name = "${assembly.asm_name}".replace(" ", "_")
  def asm_last = "${assembly.ftp_path}".split('/').last()
  """
  echo LAST ${asm_last}
  echo ${assembly.ftp_path}/${asm_last}_protein.faa.gz -o ${assembly.assembly_accession}.faa.gz
  curl  ${assembly.ftp_path}/${asm_last}_protein.faa.gz -o ${assembly.assembly_accession}.faa.gz
  """
}

process execute_interproscan {

  publishDir 'annotation/interproscan', mode: 'copy', overwrite: false

  input:
  file(seq) from assembly_faa

  output:
  file '*gff3' into interpro_gff3
  file '*html.tar.gz' into interpro_html
  file '*svg.tar.gz' into interpro_svg
  file '*tsv' into interpro_tsv
  file '*xml' into interpro_xml

  script:
  n = seq.name
  """
  interproscan.sh -appl ProDom,HAMAP,TIGRFAM,SUPERFAMILY,PRINTS,PIRSF,COILS,ProSiteProfiles,PfamA,SMART,Phobius,SMART,SignalP_GRAM_NEGATIVE -goterms -pa -f TSV -i ${n} -d . -T . -iprlookup
  """
}
