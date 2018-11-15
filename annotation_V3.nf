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
log.info "====================================="
log.info "input                  : ${params.input}"
log.info "Blast COG              : ${params.blast_cog}"
log.info "Orthofinder            : ${params.orthofinder}"
log.info "Orthofinder path       : ${params.genome_faa_folder}"


  Channel
    .fromPath(params.input)
    .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
    .set { faa_genomes }


process prepare_orthofinder {

  echo true

  conda 'bioconda::orthofinder=2.2.7'

  input:
    file genome_list from faa_genomes.collect()

  output:
    file 'Results_Nov*/WorkingDirectory/Species*.fa' into species_fasta
    file 'Results_Nov*/WorkingDirectory/BlastDBSpecies*.phr' into species_blastdb
    file 'Result*' into result_dir

  script:
  """
  orthofinder -op -a 8 -f . > of_prep.tab
  """
}

process blast {

  echo true

  input:
  file complete_dir from result_dir
  each seq from species_fasta
  each blastdb from species_blastdb

  output:
  file "${complete_dir.baseName}/Blast_${species_1}_${species_2}.txt" into blast_results

  script:
  blastdb_name = blastdb.getBaseName()
  blastdb_path = blastdb.getParent()
  seq_name = seq.getBaseName()
  species_1 =  (seq_name =~ /Species(\d+)/)[0][1]
  species_2 =  (blastdb_name =~ /BlastDBSpecies(\d+)/)[0][1]

  """
  blastp -outfmt 6 -evalue 0.001 -query $seq -db $blastdb_path/$blastdb_name > ${complete_dir}/WorkingDirectory/Blast_${species_1}_${species_2}.txt
  """

}

blast_results.collect().set { blast_files }

process print_cmds {
    echo true

    input:
    val blast_cmd from blast_files

    script:
    """
    echo "aaa ${blast_cmd} bbb"
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
