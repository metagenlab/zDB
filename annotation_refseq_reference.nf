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

  conda 'bioconda::biopython=1.68'

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
no_mapping_uniparc_records = []

lines = []
for record in records:

    checksum = CheckSum.seguid(record.seq)
    print(checksum)
    sql1 = 'select uniparc_accession from uniparc.uniparc_accession where sequence_hash=?'
    cursor.execute(sql1, (checksum,))
    hits = cursor.fetchall()

    # check if uniparc match

    if len(hits) == 0:
        no_mapping_uniparc_records.append(record)
    else:
        # fetch interpro annotation
        sql = 'select uniparc_accession,t0.uniparc_id,entry_accession,entry_name,accession,description,evidence_name,db_name,name,start,end,score from uniparc.uniparc_accession t0 inner join uniparc_match t1 on t0.uniparc_id=t1.uniparc_id inner join match_entry t2 on t1.entry_id=t2.entry_id inner join interpro_entry t3 on t1.interpro_id=t3.interpro_id inner join evidence t4 on t1.evidence_id=t4.evidence_id inner join database t5 on t2.db_id=t5.db_id inner join interpro_type t6 on t3.type_id=t6.type_id where t0.sequence_hash=? and t5.db_name!="MOBIDBLT"'
        cursor.execute(sql, (checksum,))
        hits = cursor.fetchall()

        # some uniparc entries wont have any interpro annotation 
        # ATTENTION: should be the same uniparc release: if uniparc release is more recent than the 
        # annotated uniparc release from InterPro, new uniparc entries could be missing from the 
        # interpro annotation
        if len(hits) != 0:
          for hit in hits:
            lines.append("%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t\\n" % (record.id,
                                                                                                        hit[0],
                                                                                                        hit[1],
                                                                                                        hit[2],
                                                                                                        hit[3],
                                                                                                        hit[4],
                                                                                                        hit[5],
                                                                                                        hit[6],
                                                                                                        hit[7],
                                                                                                        hit[8],
                                                                                                        hit[9],
                                                                                                        hit[10],
                                                                                                        hit[11]))
uniparc_map.writelines(lines)

SeqIO.write(no_mapping_uniparc_records, no_uniparc_mapping, "fasta")
  """
}


no_uniparc_mapping_faa.collectFile(name: 'merged.faa', newLine: true)
    .into { merged_no_uniparc_faa }


process execute_interproscan_no_uniparc_matches {

  publishDir 'refseq_annotation/interproscan', mode: 'link', overwrite: true

  cpus 16
  memory '16 GB'
  container 'metagenlab/annotation-pipeline:1.2.1'

  input:
  file(seq) from merged_no_uniparc_faa.splitFasta( by: 5000, file: "no_uniparc_match_chunk_" )

  output:
  file '*tsv' into interpro_tsv_no_uniparc

  script:
  n = seq.name
  """
  echo $INTERPRO_HOME/interproscan.sh -f TSV -i ${n} -d . -T . --disable-precalc -cpu ${task.cpus}
  bash $INTERPRO_HOME/interproscan.sh -f TSV -i ${n} -d . -T . --disable-precalc -cpu ${task.cpus} >> ${n}.log
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