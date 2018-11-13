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
params.databases_dir = "$PWD/databases"
log.info "====================================="
log.info "input                  : ${params.input}"

Channel
  .fromPath(params.input)
  .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
  .splitFasta( by: 1000 )
  .set { faa_chunks }

  process blastThemAll {

      input:
      file 'seq' from faa_chunks
      output:
      file 'blast_result'
      script:
      """
      blastp -db $params.databases_dir/COG/prot2003-2014.fa -query seq -outfmt 6 > blast_result
      """
  }

  blast_result
      .collectFile(name: "$PWD/blast_COG/table.tab")
