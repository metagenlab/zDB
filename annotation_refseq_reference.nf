#!/usr/bin/env nextflow

log.info params.input
params.databases_dir = "$PWD/databases"
params.kofamscan = false
params.mapping_uniparc = true

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
                   faa_list_2}


process get_uniparc_mapping {

  conda 'bioconda::biopython=1.68'

  publishDir 'refseq_annotation/uniparc_mapping/', mode: 'link', overwrite: true
  cpus 1

  echo = true

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

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum
import glob
import re

conn = sqlite3.connect("${params.databases_dir}/interproscan/annotated_uniparc/uniparc_match.db")
cursor = conn.cursor()

sql = 'attach "${params.databases_dir}/uniprot/uniparc/uniparc.db" as uniparc'
cursor.execute(sql,)

fasta_file = glob.glob("*.faa")[0]

basename = re.sub(".faa", "", fasta_file)

uniparc_map = open('%s_uniparc_mapping_interpro.tab' % basename, 'w')
no_uniparc_mapping = open('%s_no_uniparc_mapping.faa' % basename, 'a')

uniparc_map.write("locus_tag\\tuniparc_accession\\tuniparc_id\\tentry_accession\\tentry_name\\taccession\\tdescription\\tevidence_name\\tdb_name\\tname\\tstart\\tend\\tscore\\n")

records = SeqIO.parse(fasta_file, "fasta")
no_mapping_uniprot_records = []
no_mapping_uniparc_records = []

for record in records:
    sql = 'select uniparc_accession,t0.uniparc_id,entry_accession,entry_name,accession,description,evidence_name,db_name,name,start,end,score from uniparc.uniparc_accession t0 inner join uniparc_match t1 on t0.uniparc_id=t1.uniparc_id inner join match_entry t2 on t1.entry_id=t2.entry_id inner join interpro_entry t3 on t1.interpro_id=t3.interpro_id inner join evidence t4 on t1.evidence_id=t4.evidence_id inner join database t5 on t2.db_id=t5.db_id inner join interpro_type t6 on t3.type_id=t6.type_id where t0.sequence_hash=?'
    cursor.execute(sql, (CheckSum.seguid(record.seq),))
    hits = cursor.fetchall()
    if len(hits) == 0:
        no_mapping_uniparc_records.append(record)
        no_mapping_uniprot_records.append(record)
    else:
        uniparc_map.write("%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t\\n" % (record.id,
                                                                                             hits[0][0],
                                                                                             hits[0][1],
                                                                                             hits[0][2],
                                                                                             hits[0][3],
                                                                                             hits[0][4],
                                                                                             hits[0][5],
                                                                                             hits[0][6],
                                                                                             hits[0][7],
                                                                                             hits[0][8],
                                                                                             hits[0][9],
                                                                                             hits[0][10],
                                                                                             hits[0][11]))

SeqIO.write(no_mapping_uniparc_records, no_uniparc_mapping, "fasta")
  """
}

no_uniparc_mapping_faa.collectFile(name: 'merged.faa', newLine: true)
    .into { merged_no_uniparc_faa }

process execute_interproscan_no_uniparc_matches {

  publishDir 'refseq_annotation/interproscan', mode: 'link', overwrite: true

  cpus 8
  memory '16 GB'
  conda 'anaconda::openjdk=8.0.152'

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

process execute_kofamscan {

  publishDir 'refseq_annotation/KO', mode: 'link', overwrite: true

  cpus 4
  memory '8 GB'

  when:
  params.kofamscan == true

  input:
  file "*" from faa_list_2

  output:
  file '*tab'

  script:
  """
  export PATH="$PATH:/home/tpillone/work/dev/annotation_pipeline_nextflow/bin/KofamScan/"
  input_file=`ls *faa`
  echo \$input_file
  echo exec_annotation \${input_file} -p ${params.databases_dir}/kegg/profiles/prokaryote.hal -k ${params.databases_dir}/kegg/ko_list --cpu ${task.cpus} -o \${input_file/faa/tab}
  exec_annotation \${input_file} -p ${params.databases_dir}/kegg/profiles/prokaryote.hal -k ${params.databases_dir}/kegg/ko_list --cpu ${task.cpus} -o \${input_file/faa/tab}
  """
}
