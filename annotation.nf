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

log.info "====================================="
log.info "input                  : ${params.input}"

Channel
  .fromPath(params.input)
  .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
  .set { datasets }
