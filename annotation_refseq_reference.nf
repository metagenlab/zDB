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
  sleep 1
  """
}


process extract_gz {

  publishDir 'refseq_annotation/refseq_faa', mode: 'link', overwrite: false

  input:
    file assembly_faa_gz

  output:
    file "${assembly_faa_gz.baseName}" into assembly_faa

  script:
  """
  zcat ${assembly_faa_gz} > ${assembly_faa_gz.baseName}
  """
}

// sort uniparc sequences from non-uniparc sequences

//assembly_faa.flatten().map { it }.filter { (it.text =~ /(>)/).size() < 1000 }.set { small_genomes }

assembly_faa.flatten().into {faa_list_1
                            faa_list_2}

/*
assembly_faa.collectFile(name: 'merged.faa', newLine: true)
    .into { faa_list_1
            faa_list_2 }
*/

process get_uniparc_mapping {

  conda 'bioconda::biopython=1.68'

  publishDir 'refseq_annotation/uniparc_mapping/', mode: 'link', overwrite: true

  when:
  params.uniparc == true

  input:
  file "${x}.faa" from faa_list_1

  output:
  file "${x}_uniparc_mapping.tab" into uniparc_mapping_tab
  file "${x}_no_uniparc_mapping.faa" into no_uniparc_mapping_faa
  file "${x}_uniparc_mapping.faa" into uniparc_mapping_faa

  script:

  """
#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum

conn = sqlite3.connect("${params.databases_dir}/uniprot/uniparc/uniparc.db")
cursor = conn.cursor()

fasta_file = "${x}.faa"

uniparc_map = open('uniparc_mapping.tab', 'w')
no_uniparc_mapping = open('no_uniparc_mapping.faa', 'a')
uniparc_mapping_faa = open('uniparc_mapping.faa', 'a')

uniparc_map.write("locus_tag\\tuniparc_id\\tuniparc_accession\\tstatus\\n")
uniprot_map.write("locus_tag\\tuniprot_accession\\ttaxon_id\\tdescription\\n")

records = SeqIO.parse(fasta_file, "fasta")
no_mapping_uniprot_records = []
no_mapping_uniparc_records = []
mapping_uniparc_records = []

for record in records:
    sql = 'select uniparc_id,uniparc_accession where sequence_hash=?'
    cursor.execute(sql, (CheckSum.seguid(record.seq),))
    hits = cursor.fetchall()
    if len(hits) == 0:
        no_mapping_uniparc_records.append(record)
        no_mapping_uniprot_records.append(record)
    else:
        mapping_uniparc_records.append(record)
        uniparc_map.write("%s\\t%s\\t%s\\t%s\\n" % (record.id,
                                                    hits[0][0],
                                                    hits[0][1]))

SeqIO.write(no_mapping_uniparc_records, no_uniparc_mapping, "fasta")
SeqIO.write(mapping_uniparc_records, uniparc_mapping_faa, "fasta")
  """
}

no_uniparc_mapping_faa.collectFile(name: 'merged.faa', newLine: true)
    .into { merged_no_uniparc_faa }

process execute_interproscan_no_uniparc_matches {

  publishDir 'refseq_annotation/interproscan', mode: 'link', overwrite: true

  cpus 8
  memory '16 GB'
  conda 'anaconda::openjdk=8.0.152'

  when:
  params.interproscan == true

  input:
  file(seq) from merged_no_uniparc_faa.splitFasta( by: 500, file: "no_uniparc_match_chunk_" )

  output:
  file '*gff3' into interpro_gff3_no_uniparc
  file '*html.tar.gz' into interpro_html_no_uniparc
  file '*svg.tar.gz' into interpro_svg_no_uniparc
  file '*tsv' into interpro_tsv_no_uniparc
  file '*xml' into interpro_xml_no_uniparc
  file '*log' into interpro_log_no_uniparc

  script:
  n = seq.name
  """
  echo $INTERPRO_HOME/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . --disable-precalc -cpu ${task.cpus}
  bash $INTERPRO_HOME/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . --disable-precalc -cpu ${task.cpus} >> ${n}.log
  """
}

/*
process interproscan_annotation_uniparc_matches {

  publishDir 'refseq_annotation/interproscan', mode: 'link', overwrite: true

  cpus 8
  memory '8 GB'
  conda 'anaconda::openjdk=8.0.152'

  when:
  params.interproscan == true

  input:
  file(seq) from uniparc_mapping_faa.splitFasta( by: 1000, file: "uniparc_match_chunk_" )

  output:
  file '*gff3' into interpro_gff3_uniparc
  file '*html.tar.gz' into interpro_html_uniparc
  file '*svg.tar.gz' into interpro_svg_uniparc
  file '*tsv' into interpro_tsv_uniparc
  file '*xml' into interpro_xml_uniparc
  file '*log' into interpro_log_uniparc

  script:
  n = seq.name
  """
  echo $INTERPRO_HOME/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . -iprlookup -cpu ${task.cpus}
  bash $INTERPRO_HOME/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . -iprlookup -cpu ${task.cpus} >> ${n}.log
  """
}
*/

process execute_kofamscan {

  publishDir 'refseq_annotation/KO', mode: 'link', overwrite: true

  cpus 4
  memory '8 GB'

  when:
  params.ko == true

  input:
  file "${x}.faa" from faa_list_2

  output:
  file '*tab'

  script:
  n = seq.name
  """
  export PATH="$PATH:/home/tpillone/work/dev/annotation_pipeline_nextflow/bin/KofamScan/"
  exec_annotation ${x}.faa -p ${params.databases_dir}/kegg/profiles/prokaryote.hal -k ${params.databases_dir}/kegg/ko_list --cpu ${task.cpus} -o ${x}.tab
  """
}
