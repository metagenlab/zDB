#!/usr/bin/env nextflow
/*
 * Author:
 * - Trestan Pillonel <trestan.pillonel@gmail.com>
 *
 */


/*
 * Pipeline input params
 */

params.input = "faa/*.faa" 	// input sequences
log.info params.input
params.databases_dir = "$PWD/databases"
params.blast_cog = true
params.orthofinder = true
params.genome_faa_folder = "$PWD/faa"

params.local_sample_sheet = "local_assemblies.tab"
params.ncbi_sample_sheet = "ncbi_assemblies.tab"


log.info "====================================="
log.info "input                  : ${params.input}"
log.info "Blast COG              : ${params.blast_cog}"
log.info "Orthofinder            : ${params.orthofinder}"
log.info "Orthofinder path       : ${params.genome_faa_folder}"

//assembly_accession_list = ""

// Each Sample
if (params.ncbi_sample_sheet != false){
  Channel.fromPath( file(params.ncbi_sample_sheet) )
                      .splitCsv(header: true, sep: '\t')
                      .map{row ->
                          // get the list of accessions
                          def assembly_accession = row."Genbank"
                          return "${assembly_accession}"
                      }
                      .into{
                          assembly_accession_list
                      }
}
if (params.local_sample_sheet != false){
  Channel.fromPath( file(params.local_sample_sheet) )
                      .splitCsv(header: true, sep: '\t')
                      .map{row ->
                          // get the list of accessions
                          def gbk_path = row."gbk_path"
                          return gbk_path
                      }
                      .map { file(it) }
                      .set { local_gbk_list }
}

// only define process if nedded
if (params.local_sample_sheet != false){
  process copy_local_assemblies {

    publishDir 'data/gbk_local', mode: 'copy', overwrite: true

    cpus 1

    when:
    params.local_sample_sheet != false

    input:
    file(local_gbk) from local_gbk_list

    output:
    file "${local_gbk.name}.gz" into raw_local_gbffs

    script:

    """
    gzip -f ${local_gbk.name}
    """
  }
}

// only define process if nedded
if (params.ncbi_sample_sheet != false){
  process download_assembly {

    conda 'bioconda::biopython=1.68'

    publishDir 'data/gbk_ncbi', mode: 'copy', overwrite: true

    when:
    params.ncbi_sample_sheet != false

    input:
    val assembly_accession_list from assembly_accession_list.collect()

    cpus 1

    output:
    file '*.gbff.gz' into raw_ncbi_gbffs

    script:
    //accession = assembly_accession[0]
    """
  #!/usr/bin/env python

  import re
  from ftplib import FTP
  from Bio import Entrez, SeqIO
  Entrez.email = "trestan.pillonel@chuv.ch"
  Entrez.api_key = "719f6e482d4cdfa315f8d525843c02659408"

  accession_list = "${assembly_accession_list}".split(' ')
  for accession in accession_list:
    handle1 = Entrez.esearch(db="assembly", term="%s" % accession)
    record1 = Entrez.read(handle1)

    ncbi_id = record1['IdList'][-1]
    print(ncbi_id)
    handle_assembly = Entrez.esummary(db="assembly", id=ncbi_id)
    assembly_record = Entrez.read(handle_assembly, validate=False)
    ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
    print(ftp_path)
    ftp=FTP('ftp.ncbi.nih.gov')
    ftp.login("anonymous","trestan.pillonel@unil.ch")
    ftp.cwd(ftp_path)
    filelist=ftp.nlst()
    filelist = [i for i in filelist if 'genomic.gbff.gz' in i]
    print(filelist)
    for file in filelist:
      ftp.retrbinary("RETR "+file, open(file, "wb").write)
    """
  }
}

// merge local and ncbi gbk into a single channel
if (params.ncbi_sample_sheet != false && params.local_sample_sheet == false) {
println "ncbi"
raw_ncbi_gbffs.collect().into{all_raw_gbff}
}
else if(params.ncbi_sample_sheet == false && params.local_sample_sheet != false) {
println "local"
raw_local_gbffs.collect().into{all_raw_gbff}
}
else {
println "both"
raw_ncbi_gbffs.mix(raw_local_gbffs).into{all_raw_gbff}
}

process gbk_check {

  publishDir 'data/gbk_edited', mode: 'copy', overwrite: true

  conda 'bioconda::biopython=1.68'

  cpus 2

  input:
  file(all_gbff) from all_raw_gbff.collect()

  output:
  file "*merged.gbk" into edited_gbks

  script:
  println all_gbff
  """
  gbff_check.py -i ${all_gbff} -l 1000
  """
}

process convert_gbk_to_faa {

  publishDir 'data/faa_locus', mode: 'copy', overwrite: true

  conda 'bioconda::biopython=1.68'

  cpus 2

  input:
  each file(edited_gbk) from edited_gbks

  output:
  file "*.faa" into faa_files

  script:
  """
#!/usr/bin/env python
print("${edited_gbk}")
from Bio import Entrez, SeqIO

records = SeqIO.parse("${edited_gbk}", 'genbank')
edited_records = open("${edited_gbk.baseName}.faa", 'w')
for record in records:
  for feature in record.features:
      if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers:
          feature.name = feature.qualifiers["locus_tag"]
          edited_records.write(">%s %s\\n%s\\n" % (feature.qualifiers["locus_tag"][0],
                                                   record.description,
                                                   feature.qualifiers['translation'][0]))
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
