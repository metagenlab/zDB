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
params.blast_cog = false
params.orthofinder = true
params.genome_faa_folder = "$PWD/faa"
log.info "====================================="
log.info "input                  : ${params.input}"
log.info "Blast COG              : ${params.blast_cog}"
log.info "Orthofinder            : ${params.orthofinder}"
log.info "Orthofinder path       : ${params.genome_faa_folder}"

Channel
  .fromPath(params.input)
  .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
  .splitFasta( by: 1000 )
  .set { faa_chunks }

genome_folder = Channel
              .from(params.genome_faa_folder)
/*
process blast_COG {

  conda 'bioconda::blast=2.7.1'

  when:
  params.blast_cog == true

  input:
  file 'seq' from faa_chunks

  output:
  set file('blast_result') into merged_blasts

  script:
  """
  blastp -db $params.databases_dir/COG/prot2003-2014.fa -query seq -outfmt 6 > blast_result
  """
}
 */
process orthofinder_preparation {

  conda 'bioconda::orthofinder=2.2.7'

  input:
  val faa_folder from genome_folder

  output:
  stdout orthofinder_prep_output

  when:
  params.orthofinder == true

  script:
  log.info "-----------> ${params.genome_faa_folder}"
  myDir = file("${params.genome_faa_folder}/Results*/WorkingDirectory/*")

  if (myDir.size() > 0){

    """
    rm -rf  ${params.genome_faa_folder}/Results*
    orthofinder -op -a 8 -f ${params.genome_faa_folder} > orthofinder_prep_output
    """
  }
  else{

  """
  orthofinder -op -a 8 -f ${params.genome_faa_folder}
  """
}
}

process write_data {

    cpus 2

    memory '3 GB'

    input:
    file orthofinder from orthofinder_prep_output

"""
#!/usr/bin/env python

with open("$orthofinder", 'r') as f:
  for n, row in enumerate(f):
    print (n, row)
"""
}
