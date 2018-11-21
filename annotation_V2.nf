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
  .splitFasta( by: 1000 )
  .set { faa_chunks }

genome_folder = Channel
              .from(params.genome_faa_folder)

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


process orthofinder_preparation {

  conda 'bioconda::orthofinder=2.2.7'

  input:
  val faa_folder from genome_folder

  output:
  file 'of_prep.tab' into orthofinder_prep_output

  when:
  params.orthofinder == true

  script:
  myDir = file("${params.genome_faa_folder}/Results*/WorkingDirectory/*")

  if (myDir.size() > 0){
    log.info "removing existing dir..."
    """
    rm -rf  ${params.genome_faa_folder}/Results*
    orthofinder -op -a 8 -f ${params.genome_faa_folder} > of_prep.tab
    """
  }
  else{
    log.info "first attempt..."
    """
    orthofinder -op -a 8 -f ${params.genome_faa_folder} > of_prep.tab
    """
  }
}

/*
process parse_orthofinder_blast {
    echo true

    input:
    file 'of_prep.tab' from orthofinder_prep_output2

    output:
    val blast_cmd

    println 'executing script'
    """
    #!/usr/bin/env python3
    import sys
    print('ok')
    x = 0
    for line in open("of_prep.tab", 'r'):
        if 'blastp -outfmt ' in line:
          blast_cmd = line.rstrip()
          x+=1
    print(x)
    """


    for( line in x.readLines() ) {
        if ('blast' in line){println line}
    }


}

.splitText()
                       .subscribe { print it }


*/


orthofinder_prep_output.collectFile().splitText().set { blast_cmds }

process print_cmds {
    echo true

    input:
    val blast_cmd from blast_cmds

    println "line:"
    println blast_cmd.class

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
