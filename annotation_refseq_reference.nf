#!/usr/bin/env nextflow

log.info "====================================="
log.info "Database folder        : ${params.databases_dir}"

def get_interpro_version(){
    def url = 'https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload';
    def page = url.toURL().text;
    def version = (page =~ /interproscan-(\d.\d+-\d+.\d)/);
    return version[0][0];
}


process download_refseq_table {

  publishDir 'refseq_annotation/genome_table/', mode: 'link', overwrite: false

  output:
    file 'assembly_summary_refseq.txt' into refseq_table
    file 'refseq_release.txt' into refseq_release

  script:
  """
  curl -S ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt | tail -n +2 | sed 's/# assembly_accession/assembly_accession/' > assembly_summary_refseq.txt
  curl ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER -o refseq_release.txt
  """
}

/* .flatten().map { it }.filter { (it.refseq_category ==~ "reference genome")} */
refseq_table.splitCsv(sep: '\t', header: true).filter { (it.refseq_category ==~ "reference genome|representative genome")}.set { representative_and_reference_genomes }

process download_reference_genomes {

  maxForks 2
  maxRetries 3
  errorStrategy 'ignore'

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
  sleep 0.5
  """
}


process extract_gz {

  publishDir 'refseq_annotation/refseq_faa', mode: 'link', overwrite: false
  maxForks 1

  input:
    file assembly_faa_gz

  output:
    file "${assembly_faa_gz.baseName}" into assembly_faa

  script:
  """
  zcat ${assembly_faa_gz} > ${assembly_faa_gz.baseName}
  """
}



assembly_faa.into {faa_list_1
                   faa_list_2
                   faa_list_3}


process get_uniparc_mapping {

  //conda 'bioconda::biopython=1.68 anaconda::pandas'

  publishDir 'refseq_annotation/uniparc_mapping/', mode: 'link', overwrite: true
  cpus 1

  when:
  params.mapping_uniparc == true

  input:
  file "*" from faa_list_1

  output:
  file "*_uniparc_mapping_interpro.tab" into uniparc_mapping_tab
  file "*_no_uniparc_mapping.faa" into no_uniparc_mapping_faa

  script:

  """
#!/usr/bin/env python
import refseq 
import glob
print("ok4")
fasta_file = glob.glob("*.faa")[0]
refseq.uniparc_interpro_annotations(fasta_file, "${params.databases_dir}")
  """
}


no_uniparc_mapping_faa.collectFile(name: 'merged.faa', newLine: true)
    .set { merged_no_uniparc_faa }


process execute_interproscan_no_uniparc_matches {

  publishDir 'refseq_annotation/interproscan', mode: 'link', overwrite: true

  cpus 16
  memory '16 GB'
  container "$params.container_interproscan"

  input:
  file(seq) from merged_no_uniparc_faa.splitFasta( by: 5000, file: "no_uniparc_match_chunk_" )

  output:
  file '*tsv' into interpro_tsv_no_uniparc

  script:
  n = seq.name
  """
  echo ${params.interproscan_home}/interproscan.sh -f TSV -i ${n} -d . -T . --disable-precalc -cpu ${task.cpus}
  bash ${params.interproscan_home}/interproscan.sh -f TSV -i ${n} -d . -T . --disable-precalc -cpu ${task.cpus} --applications TIGRFAM,SFLD,SMART,CDD,SUPERFAMILY,PANTHER,Gene3D,PIRSF,Pfam >> ${n}.log
  """
}

process execute_kofamscan {

  publishDir 'refseq_annotation/KO', mode: 'link', overwrite: true

  container 'quay.io/biocontainers/kofamscan:1.3.0--0'

  cpus 4
  memory '8 GB'

  input:
  file "*" from faa_list_2

  output:
  file '*tab'

  script:
  """
  input_file=`ls *faa`
  echo exec_annotation \${input_file} -p ${params.databases_dir}/kegg/profiles/prokaryote.hal -k ${params.databases_dir}/kegg/ko_list --cpu ${task.cpus} -o \${input_file/faa/tab}
  exec_annotation \${input_file} -p ${params.databases_dir}/kegg/profiles/prokaryote.hal -k ${params.databases_dir}/kegg/ko_list --cpu ${task.cpus} -o \${input_file/faa/tab}
  """
}


process execute_rpsblast_COG {

  publishDir 'refseq_annotation/COG', mode: 'link', overwrite: true

  container 'quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4'

  cpus 4
  memory '2 GB'

  input:
  file (genome) from faa_list_3

  output:
  file "${genome.baseName}.tab" into blast_result

  script:
  """
  rpsblast -db $params.databases_dir/cdd/profiles/Cog -query ${genome} -outfmt 6 -evalue 0.001 -num_threads ${task.cpus} > ${genome.baseName}.tab
  """
}